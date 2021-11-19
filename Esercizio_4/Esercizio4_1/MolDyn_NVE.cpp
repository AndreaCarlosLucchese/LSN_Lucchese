/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <iomanip>
#include "MolDyn_NVE.h"

using namespace std;

int main(){ 

  cout << "If you want to use old configuation type 1 " << endl;
  cin >> appo2;
  cout << "Do you want to rescale the coordinates in order to match the desired T?" << endl;
  cout << "Type 1 to confirm, 0 to decline:" << endl;
  cin >> appo;
  cout << "How many blocks ?" << endl;
  cin >> nblk;
  Input();             //Inizialization
  int nconf = 1;
  for(int iblk=1;iblk<=nblk;iblk++){
     reset(iblk);
     cout << "Blocco numero " << iblk << endl;
      for(int istep=0;istep<nstep/10;istep++){
      Move();
      Measure();
      accumulate();
    }
    Average(iblk);
  }
    ConfFinal();  
  //Calcolo il g(r) medio
  return 0;
  
}


void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;
  double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open("input.dat"); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1;
  it = 2;
  ie = 3;
  n_props = 4; //Number of observables

//Read initial configuration
  cout << "Read initial configuration from file config.0 " << endl << endl;
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();

//Prepare initial velocities
   cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
   double sumv[3] = {0.0, 0.0, 0.0};
   for (int i=0; i<npart; ++i){
     vx[i] = rand()/double(RAND_MAX) - 0.5;
     vy[i] = rand()/double(RAND_MAX) - 0.5;
     vz[i] = rand()/double(RAND_MAX) - 0.5;

     sumv[0] += vx[i];
     sumv[1] += vy[i];
     sumv[2] += vz[i];
   }
   for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
   double sumv2 = 0.0, fs;
   for (int i=0; i<npart; ++i){
     vx[i] = vx[i] - sumv[0];
     vy[i] = vy[i] - sumv[1];
     vz[i] = vz[i] - sumv[2];

     sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
   }
   sumv2 /= (double)npart;

   if (appo2==0) {
   fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
   for (int i=0; i<npart; ++i){
     vx[i] *= fs;
     vy[i] *= fs;
     vz[i] *= fs;

     xold[i] = Pbc(x[i] - vx[i] * delta);
     yold[i] = Pbc(y[i] - vy[i] * delta);
     zold[i] = Pbc(z[i] - vz[i] * delta);
   }
   
  }

  else{
    ReadConf.open("old.0");
    for (int i=0; i<npart; ++i){
    ReadConf >> xold[i];
    ReadConf >> yold[i];
    ReadConf >> zold[i];
   }
    
  }

  ReadConf.close();

    if (appo == 1){
      double v2tot = 0;
    for (int i=0; i<npart; ++i){
        vx[i] = (x[i] - xold[i])/delta;
        vy[i] = (y[i] - yold[i])/delta;
        vz[i] = (z[i] - zold[i])/delta;

        v2tot += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]; //trovo il v^2 di ogni particella
      }
      double v2tot_new = 3*temp*npart;
      double cf = v2tot_new / v2tot; //trovo il fattore di correzione

    for (int i=0; i<npart; ++i){ //adesso riscalo
        xold[i] = Pbc(x[i] - sqrt(cf)* vx[i] *delta);
        yold[i] = Pbc(y[i] - sqrt(cf)* vy[i] *delta);
        zold[i] = Pbc(z[i] - sqrt(cf)* vz[i] *delta);

    }
  } 
   return;
}

void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}

void Measure(){ //Properties measurement
  double v, tij, k, vij;
  double dx, dy, dz, dr;
  ofstream epot,ekin,temp,etot;
  epot.open("output_epot.dat",ios::app);
  ekin.open("output_ekin.dat",ios::app);
  etot.open("output_etot.dat",ios::app);
  temp.open("output_temp.dat",ios::app);
  v = 0.0; //reset observables
  k = 0.0;
//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
     dy = Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
     dz = Pbc( zold[i] - zold[j] ); // => EPot should be computed with r(t)

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);
    //update of the histogram of g(r)

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
       v += vij;
    }          
  }
  }
//Kinetic energy
    for (int i=0; i<npart; ++i) tij += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
    k += tij;

    walker_av[iv] = v;  // /(double)npart; //Potential energy per particle
    walker_av[ik] = k;// /(double)npart;
    walker_av[it] = (2.0/3.0)*k;
    walker_av[ie] = k+v;
  //  Epot << walker_av[iv]  << endl;
  // Pres << walker_av[iw] << endl;
    epot << v/(double)npart << endl;
    etot << (k+v)/(double)npart << endl;
    ekin << k/(double)npart<< endl;
    temp << (2.0/3.0)*k/(double)npart << endl;

    epot.close();
    ekin.close();
    temp.close();
    etot.close();


    return;
}

void accumulate(){
  for(int i=0; i<4; i++){
    blk_av[i]=blk_av[i]+walker_av[i];
  }
  norm=norm+1;
}

void reset(int iblk){
  for(unsigned int i=0; i<4 ; i++){
    blk_av[i]=0;
  }
  norm=0;
}

void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();

  WriteConf.open("old.final");

  for (int i=0; i<npart; ++i){
    WriteConf << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
  }
  WriteConf.close();


  return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}

double Error(double sum, double sum2, int iblk)
{
    if( iblk == 1 ) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

void Average(int iblk){

  double err_pot, err_k, err_temp, err_etot;
  ofstream epot,ekin,temp,etot;
  epot.open("ave_epot.out",ios::app);
  ekin.open("ave_ekin.out",ios::app);
  etot.open("ave_etot.out",ios::app);
  temp.open("ave_temp.out",ios::app);

  const int wd=12;

    stima_pot = blk_av[iv]/norm/(double)npart;  //Potential energy
    global_m[iv] += stima_pot;
    global_m2[iv] += stima_pot*stima_pot;
    err_pot=Error(global_m[iv],global_m2[iv],iblk);

    stima_ekin = blk_av[ik]/norm/(double)npart;  //Potential energy
    global_m[ik] += stima_ekin;
    global_m2[ik] += stima_ekin*stima_ekin;
    err_k=Error(global_m[ik],global_m2[ik],iblk); 

    stima_temp = blk_av[it]/norm/(double)npart;  //Potential energy
    global_m[it] += stima_temp;
    global_m2[it] += stima_temp*stima_temp;
    err_temp=Error(global_m[it],global_m2[it],iblk);
    
    stima_etot = blk_av[ie]/norm/(double)npart;  //Potential energy
    global_m[ie] += stima_etot;
    global_m2[ie] += stima_etot*stima_etot;
    err_etot=Error(global_m[ie],global_m2[ie],iblk);


    //Potential energy per particle
    epot << setw(wd) << iblk <<  setw(wd) << stima_pot << setw(wd) << global_m[iv]/(double)iblk << setw(wd) << err_pot << endl;
    //Kinetic energy per particle
    ekin << setw(wd) << iblk <<  setw(wd) << stima_ekin << setw(wd) << global_m[ik]/(double)iblk << setw(wd) << err_k << endl;
    //Temperature per particle
    temp << setw(wd) << iblk <<  setw(wd) << stima_temp << setw(wd) << global_m[it]/(double)iblk << setw(wd) << err_temp << endl;
    //Total energy per particle
    etot << setw(wd) << iblk <<  setw(wd) << stima_etot << setw(wd) << global_m[ie]/(double)iblk << setw(wd) << err_etot << endl;

    epot.close();
    ekin.close();
    temp.close();
    etot.close();

  }



/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

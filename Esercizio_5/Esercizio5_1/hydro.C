#include "hydro.h"



Hydro::Hydro(double x, double y, double z, double step, double p){

        m_x.push_back(x);
        m_y.push_back(y);
        m_z.push_back(z);
        m_step = step;
        in_rand(m_rand);
        m_p=p;
        m_counter=0;
}

Hydro::~Hydro(){

}


double Hydro::get_radius(int i){

    return sqrt(m_x[i]*m_x[i] + m_y[i]*m_y[i] + m_z[i]*m_z[i]);
}


vector <double> Hydro::Equilibration(int M){
            for(int i=0;i<M;i++){
                metropolis_unif(i);
                m_r.push_back(get_radius(i));
            }  
    return m_r;   
}

void Hydro::reset_start(double x, double y, double z,double step) {
    m_x.clear();
    m_y.clear();
    m_z.clear();

    m_x.push_back(x);
    m_y.push_back(y);
    m_z.push_back(z);
    m_step=step;
    m_counter=0;
   // in_rand(m_rand);
}

double Hydro::prob_100(double x, double y, double z){ //funzione d'onda ground state
    double radius=sqrt(x*x + y*y+ z*z);
    return exp(-2.*radius)/M_PI;
}

double Hydro::prob_210(double x, double y, double z){ //funzione d'onda primo stato eccitato
    double radius=sqrt(x*x + y*y + z*z);
    return exp(-radius)*x*x/(32.*M_PI);
}

void Hydro::metropolis_unif(int i){ //algoritmo metropolis con distribuzione uniforme
    double x_new,y_new,z_new,a;
    if(m_p==0){
        x_new=m_x[i] +(m_rand.Rannyu()-0.5)*m_step;
        y_new=m_y[i] +(m_rand.Rannyu()-0.5)*m_step;
        z_new=m_z[i] +(m_rand.Rannyu()-0.5)*m_step;

        a=min(1.,prob_100(x_new,y_new,z_new)/prob_100(m_x[i] ,m_y[i] ,m_z[i] ));
        if(m_rand.Rannyu()<a){
            m_x.push_back(x_new);
            m_y.push_back(y_new);
            m_z.push_back(z_new);
            m_counter++;
        }
        else{
            m_x.push_back(m_x[i]);
            m_y.push_back(m_y[i]);
            m_z.push_back(m_z[i]);
        }
    }
    else{
        x_new=m_x[i] +(m_rand.Rannyu()-0.5)*m_step;
        y_new=m_y[i] +(m_rand.Rannyu()-0.5)*m_step;
        z_new=m_z[i] +(m_rand.Rannyu()-0.5)*m_step;

        a=min(1.,prob_210(x_new,y_new,z_new)/prob_210(m_x[i],m_y[i],m_z[i]));
        if(m_rand.Rannyu()<a){
            m_x.push_back(x_new);
            m_y.push_back(y_new);
            m_z.push_back(z_new);
            m_counter++;
        }
        else{
            m_x.push_back(m_x[i]);
            m_y.push_back(m_y[i]);
            m_z.push_back(m_z[i]);
        }
    }

}


void Hydro::metropolis_gauss(int i){ //algoritmo metropolis con distribuzione gaussiana
    double x_new,y_new,z_new,a;
    if(m_p==0){
        x_new=m_x[i]+m_rand.Gauss(0,1)*m_step;
        y_new=m_y[i]+m_rand.Gauss(0,1)*m_step;
        z_new=m_z[i]+m_rand.Gauss(0,1)*m_step;

        a=min(1.,prob_100(x_new,y_new,z_new)/prob_100(m_x[i],m_y[i],m_z[i]));
        if(m_rand.Rannyu()<a){
            m_x.push_back(x_new);
            m_y.push_back(y_new);
            m_z.push_back(z_new);
            m_counter++;
        }
        else{
            m_x.push_back(m_x[i]);
            m_y.push_back(m_y[i]);
            m_z.push_back(m_z[i]);
        }

    }
    else{

        x_new=m_x[i]+m_rand.Gauss(0,1)*m_step;
        y_new=m_y[i]+m_rand.Gauss(0,1)*m_step;
        z_new=m_z[i]+m_rand.Gauss(0,1)*m_step;

        a=std::min(1.,prob_210(x_new,y_new,z_new)/prob_210(m_x[i],m_y[i],m_z[i]));
        if(m_rand.Rannyu()<a){
            m_x.push_back(x_new);
            m_y.push_back(y_new);
            m_z.push_back(z_new);
            m_counter++;
        }
        else{
            m_x.push_back(m_x[i]);
            m_y.push_back(m_y[i]);
            m_z.push_back(m_z[i]);
        }

    }
}


int Hydro::get_counter(){
    return m_counter;
}

double Hydro::get_step(){  
    return m_step;
}


double* Hydro::get_pos(){
    int n=m_x.size();
    double * a= new double [3];
    a[0]=m_x[n-1];
    a[1]=m_y[n-1];
    a[2]=m_z[n-1];
    return a;
}








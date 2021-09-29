#include "lib.h"
using namespace std;

//iniziallizare classe random
void in_rand(Random &rnd){
	int seed[4];
   int p1, p2;
   ifstream Primes("../../Librerie/Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("../../Librerie/seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;
	  rnd.SaveSeed();
	  return;
}

// dist
double dist(double x1, double x2, double y1, double y2){

	return sqrt(pow(x1-x2,2)+pow(y1-y2,2));

}

void block_error(double* error, double* sum_prog, double* x, double* x2, int n){
   double* sum_prog2 = new double[n];
   
   for (int i=0; i<n; i++){
      sum_prog[i]=0;
		for (int j=0; j<i+1; j++){
         
			sum_prog[i] += x[j];
         sum_prog2[i] +=x2[j];

		}
		
		sum_prog[i] /= (i+1);
      sum_prog2[i] /= (i+1);
	}

   for(int i=0; i<n; i++){
      error[i]= sqrt( (sum_prog2[i]-pow(sum_prog[i],2)) / i);

   }


}
// Random Walk discreto

void RWD(double* RW, double* RW_error, Random &rnd,int f){
   int M=10000;
   int N=100;
   int L=M/N;

    double* err = new double [N]();
    double* sum_prog=new double [N]();
    double* ave= new double[N]();
    double* ave2= new double[N]();
    double a = 1;
    
      for(int i=0; i<N;i++){
         double appo=0;
         for(int b=0;b<L;b++){ 
            double x=0;
            double y=0;
            double z=0;
            for(int j=0;j<f;j++){
               double k= rnd.Rannyu();
               double h= rnd.Rannyu()*3;
               if(k<0.5){
                  if(h<1 and h>=0){
                     x=x-a;
                  }
                  if(h<2 and h>=1){
                     y=y-a;
                  }
                  if(h<3 and h>=2){
                     z=z-a;
                  } 
               }
               if(k>=0.5){
                  if(h<1 and h>=0){
                     x=x+a;
                  }
                  if(h<2 and h>=1){
                     y=y+a;
                  }
                  if(h<3 and h>=2){
                      z=z+a;
                  }
               }  
            }
             appo += (pow(x,2)+ pow(y,2)+ pow(z,2));
         }
      ave[i]=appo/L;
      ave2[i]=pow(ave[i],2);
   }

   block_error(err,sum_prog,ave,ave2,N);

   RW[f]=sqrt(sum_prog[N-1]);

   RW_error[f]=sqrt(err[N-1]);

   delete[] err;
   delete[] sum_prog;
   delete[] ave;
   delete[] ave2;

}

//Random Walk continuo
void RWC(double* RW, double* RW_error, Random &rnd,int f){
   int M=10000;
   int N=100;
   int L=M/N;

    double* err = new double [N]();
    double* sum_prog=new double [N]();
    double* ave= new double[N]();
    double* ave2= new double[N]();
    double a = 1;
    
      for(int i=0; i<N;i++){
         double appo=0;
         for(int b=0;b<L;b++){
            
            double x=0;
            double y=0;
            double z=0;

            for(int j=0;j<f;j++){
               double phi= rnd.theta();
               double Theta= rnd.theta()/2;
               x= x + a*sin(Theta)*cos(phi); 
               y= y + a*sin(Theta)*sin(phi);
               z= z + a*cos(Theta);
            }
            //cout << sqrt((pow(x,2)+ pow(y,2)+ pow(z,2)))<< endl;

             appo += (pow(x,2)+ pow(y,2)+ pow(z,2));

         }

      ave[i]=appo/L;
      ave2[i]=pow(ave[i],2);
   }

   block_error(err,sum_prog,ave,ave2,N);

   RW[f]=sqrt(sum_prog[N-1]);

   RW_error[f]=sqrt(err[N-1]);
   
   delete[] err;
   delete[] sum_prog;
   delete[] ave;
   delete[] ave2;


}

//Funzione utilizzate nell'esercizio 3.1 relativo all'ambito del pricing di opzioni

double d1(double S, double T, double t,double K, double r, double sigma){
      return (1/(sigma*sqrt(T-t)))*(log(S/K)+ r + pow(sigma,2)*(T-t)/2);
}

double d2(double d1, double T, double t, double sigma){
return d1 - sigma*sqrt(T-t);
}

double N(double x){
return 0.5*(1+erf(x/sqrt(2)));
}

double C(double S, double T, double t, double r, double N1, double N2, double K){

      return S*N1-K*exp(-r*(T-t))*N2;
}

double P(double S, double T, double t, double r, double N1, double N2, double K){
   return S*(N1-1)-K*exp(-r*(T-t))*(N2-1);
   }












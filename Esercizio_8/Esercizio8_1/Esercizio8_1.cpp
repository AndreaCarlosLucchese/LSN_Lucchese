#include <iostream>
#include <fstream>
#include <string>
#include "../../Librerie/random.h"
#include "../../Librerie/lib.h"
#include "Quantum.h"
using namespace std;
 
int main (int argc, char *argv[]){
   
    int M=10000;
    int N=100;
    double* x1=new double[N]();
    double* x2=new double[N]();
    double* sum_prog=new double[N]();
    double* error=new double[N]();
    double* x=new double[M]();
    int L=M/N;

    ofstream out;
    out.open("energy.txt");

    //Grid search
    // Faccio una prima griglia in cui cerco la zona in cui sono i parametri
    //Faccio una seconda griglia per migliora la stima dei parametri
    double sigma=0;
    double mu=0;
    int rep=10000;
    int d=20;
    double step=2;


    Quantum_hole test(0.,0.1,0.1,step);
    double en=test.en(test.get_radius());

    for(int i=1; i<= d;i++){       
        double t_sigma=double(1./d)*i;
        for(int j=1;j<=d;j++){
            double appo=0;
            double t_mu=double(1./d)*j;
            test.reset_start(0.,t_mu,t_sigma,step);
            for(int h=0; h< rep; h++){
                test.metropolis_unif();
                appo+=test.en(test.get_radius());
            }
            appo/=rep;
            if(appo<en){
                en=appo;
                sigma=t_sigma;
                mu=t_mu;
                cout << en << " " << mu << " " << sigma << endl;
            }
        }
    }
    double new_mu=mu;
    double new_sigma=sigma;
    for(int i=1; i<= 100;i++){       
        double t_sigma=new_sigma+0.01*i;
        for(int j=1;j<=100;j++){
            double appo=0;
            double t_mu=new_mu+0.01*j;
            test.reset_start(0.,t_mu,t_sigma,step);
            for(int h=0; h< rep; h++){
                test.metropolis_unif();
                appo+=test.en(test.get_radius());
            }
            appo/=rep;
            if(appo<en){
                en=appo;
                sigma=t_sigma;
                mu=t_mu;
            }
        }
    }

    cout << "optimal mu: " << mu << endl;
    cout<< "optimal sigma: " << sigma << endl;

    
    double new_step=5;

    ofstream dist; // calcolo l'energia e la psi
    dist.open("psi.txt");
    test.reset_start(0.,mu,sigma,new_step);

    for(int j=0; j < N; j++){
        double sum=0;
        for(int k=0; k < L; k++){
                test.metropolis_unif();
                sum+=test.en(test.get_radius());
                dist << test.get_radius() << endl;
        }
        x1[j]=sum/L;
        x2[j]=pow(x1[j],2);
    }
    block_error(error,sum_prog,x1,x2,N);
    for(int i=0; i<N;i++){
            out << sum_prog[i] << " " << error[i] << endl;
        }
    cout << "Acceptance rate:  " << (double)test.get_counter()/M << endl;
    out.close();

    dist.close();
    delete [] x1;
    delete [] x2;
    delete [] x;
    delete [] sum_prog;
    delete [] error;




    return 0;
}
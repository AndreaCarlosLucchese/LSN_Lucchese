#include <iostream>
#include <fstream>
#include <string>
#include "../../Librerie/random.h"
#include "../../Librerie/lib.h"


using namespace std;
 
int main (int argc, char *argv[]){

    Random rnd;
    in_rand(rnd);

    double S_0=100;
    double T=1;
    double K=100;
    double r=0.1;
    double sigma=0.25;
    double d1_0;
    double d2_0;
    double N1_0;
    double N2_0;
    int M=10000;
    int N1=100;
    int L=M/N1;
    double t=0;
    double S=0;
    double *C_1 = new double [N1]();
    double *P_1 = new double [N1]();
    double* C_1_2= new double[N1]();
    double* P_1_2= new double[N1]();
    double* err_C = new double [N1]();
    double* err_P= new double [N1]();
    double* sum_prog_C=new double [N1]();
    double* sum_prog2_C=new double [N1]();
    double* sum_prog_P=new double [N1]();
    double* sum_prog2_P=new double [N1]();

ofstream out;
t=T;

out.open("Econofisica.txt");
for(int i=0;i<N1;i++){
    double appo1=0;
    double appo2=0;
    for(int j=1;j<L;j++){
        S=S_0*exp((r-0.5*pow(sigma,2))*T + sigma*rnd.Gauss(0,T));
        d1_0=d1(S,T,t,K,r,sigma);
        d2_0=d2(d1_0,T,t,sigma);
        N1_0=N(d1_0);
        N2_0=N(d2_0);
        appo1 +=C(S,T,t,r,N1_0,N2_0,K);
        appo2 +=P(S,T,t,r,N1_0,N2_0,K);
    }
    C_1[i]=appo1/L;
    P_1[i]=appo2/L;
    C_1_2[i]=pow(C_1[i],2);
    P_1_2[i]=pow(P_1[i],2);

}
        block_error(err_C, sum_prog_C, C_1, C_1_2, N1);
        block_error(err_P, sum_prog_P, P_1, P_1_2, N1);

         for (int i=0; i<N1 ;i++){
         out << sum_prog_C[i] << " " << err_C[i] << " " << sum_prog_P[i] << " " << err_P[i] << endl;
        }


    out.close();

    out.open("Econofisica_2.txt");
    double dt=0.01;
    Random rnd1;
    in_rand(rnd1);
    
    for(int i=0; i < N1; i++){
        double appo1=0;
        double appo2=0;
       
        for(int k=0; k < L; k++){
            S=100;
            
            for(int j=1; j < 100; j++){
                S = S*exp((r-0.5*pow(sigma,2))*(dt)+ sigma*rnd1.Gauss(0,1)*sqrt(dt));
                }  
            
            d1_0=d1(S,T,t,K,r,sigma);
            d2_0=d2(d1_0,T,t,sigma);
            N1_0=N(d1_0);
            N2_0=N(d2_0);
            appo1 += C(S,T,t,r,N1_0,N2_0,K);
            appo2 += P(S,T,t,r,N1_0,N2_0,K);
    }
    C_1[i]=appo1/L;
    P_1[i]=appo2/L;
    C_1_2[i]=pow(C_1[i],2);
    P_1_2[i]=pow(P_1[i],2);

    }
        block_error(err_C, sum_prog_C, C_1, C_1_2, N1);
        block_error(err_P, sum_prog_P, P_1, P_1_2, N1);

         for (int i=0; i<N1 ;i++){
         out << sum_prog_C[i] << " " << err_C[i] << " " << sum_prog_P[i] << " " << err_P[i] << endl;
        }


    out.close();



delete[] C_1;
delete[] C_1_2;
delete[] P_1;
delete[] P_1_2;
delete[] err_C;
delete[] err_P;
delete[] sum_prog_C;
delete[] sum_prog_P;
delete[] sum_prog2_C;
delete[] sum_prog2_P;

return 0;
}
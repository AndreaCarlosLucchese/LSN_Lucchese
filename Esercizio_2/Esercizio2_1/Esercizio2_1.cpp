#include <iostream>
#include <fstream>
#include <string>
#include "../../Librerie/random.h"
#include "../../Librerie/lib.h"

using namespace std;
 
int main (int argc, char *argv[]){


//  Integrale con metodo della media
    Random rnd;
    in_rand(rnd);

    int M=10000;
    int N=100;
    int L=M/N;
    double* ave= new double[N]();
    double* ave2= new double[N]();
    double* err = new double [N]();
    double* sum_prog=new double [N]();
    

for(int k=0; k<N;k++){ 
    double appo =0;

    for (int i=0; i < L;i++){
        double x=0;
        x=rnd.Rannyu();
        appo += (M_PI/2)*cos(M_PI*x/2);

    }
        ave[k]= appo/L;
        ave2[k]=ave[k]*ave[k];
 }

       block_error(err,sum_prog,ave,ave2,N);

    ofstream out;
    out.open("Int_media.txt");

        for (int i=0; i<N ;i++){
         out << sum_prog[i] << " " << err[i] << endl;

        }

    out.close();


//Integrale con importance sampling
    double* ave_i= new double[N]();
    double* ave2_i= new double[N]();
    double* err_i = new double [N]();
    double* sum_prog_i=new double [N]();
    double* sum_prog2_i=new double [N]();
    
for(int k=0; k<N;k++){ 
    double appo =0;

    for (int i=0; i < L;i++){
        double x=0;
        x=1-sqrt(1-rnd.Rannyu());
        appo += (M_PI/2)*cos(M_PI*x/2)/(2-2*x);
    }
        ave_i[k]= appo/L;
        ave2_i[k]=ave_i[k]*ave_i[k];
 }

    block_error(err_i,sum_prog_i,ave_i,ave2_i,N);

    out.open("Int_samp.txt");

        for (int i=0; i<N ;i++){
         out << sum_prog_i[i] << " " << err_i[i] << endl;

        }

    out.close();

    delete [] ave;
    delete [] ave2;
    delete [] err;
    delete [] sum_prog;
    delete [] ave_i;
    delete [] ave2_i;
    delete [] err_i;
    delete [] sum_prog_i;



return 0;
}
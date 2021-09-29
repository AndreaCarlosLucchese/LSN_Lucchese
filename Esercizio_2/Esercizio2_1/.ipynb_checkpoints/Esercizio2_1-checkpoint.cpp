#include <iostream>
#include <fstream>
#include <string>
#include "lib.h"

using namespace std;
 
int main (int argc, char *argv[]){

    Random rnd;
    in_rand(rnd);

    int M=10000;
    int N=100;
    int L=M/N;
    double* ave= new double[N]();
    double* ave2= new double[N]();
    double* err = new double [N]();
    double* sum_prog=new double [N]();
    double* sum_prog2=new double [N]();
    

for(int k=0; k<N;k++){ 
    double appo =0;

    for (int i=0; i < L;i++){
        double x=0;
        x=rnd.Rannyu();
        appo += (M_PI/2)*cos(M_PI*x/2);

    }
        ave[k]= appo/L;
        ave2[k]=ave[k]*ave[k];
        cout << ave[k] << endl;
 }

        sumprog(N,sum_prog,ave);
        sumprog(N,sum_prog2,ave2);

        for (int i=0; i<N ;i++){
         err[i] = Error(sum_prog[i], sum_prog2[i], i);

        }

    ofstream out;
    out.open("Int_media.txt");

        for (int i=0; i<N ;i++){
         out << sum_prog[i] << " " << err[i] << endl;

        }

    out.close();














return 0;
}
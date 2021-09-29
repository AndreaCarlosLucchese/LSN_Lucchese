#include <iostream>
#include <fstream>
#include <string>
#include "../../Librerie/random.h"
#include "../../Librerie/lib.h"
using namespace std;
 
int main (int argc, char *argv[]){

    Random rand;
    in_rand (rand);

    double L= 0.4; //la lunghezza della corda che sto lanciando
    double d= 0.5; //la distanza tra le varie righe  
    int M=10000;
    int N=100;
    int z=50000;
    int N_hit=0;
    int L_1=M/N;
    double* ave= new double[N]();
    double* ave2= new double[N]();
    double* err = new double [N]();
    double* sum_prog=new double [N]();
//Faccio la media a blocchi del tutto, l'idea dell'algoritmo si basa sull'utilizzo di un metodo di accept reject per la generazione dell'angolo alfa
for(int k=0; k<N;k++){ 
    double appo_1=0;
    for (int i=0; i < L_1;i++){
            N_hit=0;
        for(int j=0; j<z;j++){
            double x=rand.Rannyu();
            double alfa=rand.theta();
            double appo= x+(alfa)*L;
            if (appo > x){
                if(d<=appo && d>=x){
                    N_hit=N_hit+1;
                }
            }
            if(x>=appo){
                if(d<=x && d>=appo){
                    N_hit=N_hit+1;
                }
             }
        }
        appo_1 = appo_1 + (L*z)/(N_hit*d);
    }
        ave[k]= appo_1/L_1;
        ave2[k]=ave[k]*ave[k];
 }

       
    block_error(err,sum_prog,ave,ave2,N);

    ofstream out;
    out.open("result_pi.txt");

    for (int i=0; i<N ;i++){
         out << sum_prog[i] << " " << err[i] << endl;
        }

    out.close();
    return 0;


    delete [] ave;
    delete [] ave2;
    delete [] err;
    delete [] sum_prog;


}
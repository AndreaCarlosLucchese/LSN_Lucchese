#include <iostream>
#include <fstream>
#include <string>
#include "../../Librerie/random.h"
#include "../../Librerie/lib.h"
#include "hydro.h"
using namespace std;
 
int main (int argc, char *argv[]){
    double step=0;
    int p =0;

    cout <<"Scegliere se campionare lo stato fondamentale, 0,  o il primo stato eccitato 1" << endl;

    cin >> p;

    cout << endl;
    ofstream out;
    int rep=100;
    if(p==0){
        out.open("result_100_unif.txt");
    
    }
    else{
        out.open("result_210_unif.txt");
    }
    Hydro idro(0,0,0,0,p);

   int M=10000;
    int N=100;
    double* x1=new double[N]();
    double* x2=new double[N]();
    double* sum_prog=new double[N]();
    double* error=new double[N]();
    int L=M/N;

    //Uniform Distribution


    for(int j=0; j < N; j++){
        double sum=0;
        for(int k=0; k < L; k++){
            idro.reset_start(2,1,0,6.2);
             for(int i=0;i<rep; i++){
                idro.metropolis_unif();
                
                }
                sum+=idro.get_radius();
        }
        x1[j]=sum/L;
        x2[j]=pow(x1[j],2);
    }
    block_error(error,sum_prog,x1,x2,N);
    
        for(int i=0; i<N;i++){
            out << sum_prog[i] << " " << error[i] << endl;
        }
    out.close();



// Gaussian distribution



     if(p==0){
        out.open("result_100_gauss.txt");
    
    }
    else{
        out.open("result_210_gauss.txt");
    }
    for(int j=0; j < N; j++){
        double sum=0;
        for(int k=0; k < L; k++){
            idro.reset_start(1,1,1,2.5);
             for(int i=0;i<rep; i++){
                idro.metropolis_gauss();
                
                }
                sum+=idro.get_radius();
        }
        x1[j]=sum/L;
        x2[j]=pow(x1[j],2);
    }
    block_error(error,sum_prog,x1,x2,N);
    
        for(int i=0; i<N;i++){
            out << sum_prog[i] << " " << error[i] << endl;
        }

    out.close();


   
    

    delete [] x1;
    delete [] x2;
    delete [] error;
    delete [] sum_prog;


    return 0;
}


#include <iostream>
#include <fstream>
#include <string>
#include "../../Librerie/random.h"
#include "../../Librerie/lib.h"


using namespace std;
 
int main (int argc, char *argv[]){

    Random rnd;
    in_rand(rnd);

    int N=100;

    double* RW= new double [N]();
    double* RW_error= new double [N]();
    double* RW_C= new double [N]();
    double* RW_error_C= new double [N]();
    ofstream out;


for(int j=0;j<N;j++){
    RWD(RW, RW_error,rnd,j);
}
   out.open("RW_discreto.txt");

        for (int i=0; i<N ;i++){

         out << RW[i] << " " << RW_error[i] << endl;

        }

    out.close(); 


out.open("RW_continuo.txt");
for(int j=0;j<N;j++){
    RWC(RW_C, RW_error_C,rnd,j);
}

        for (int i=0; i<N ;i++){

         out << RW[i] << " " << RW_error[i] << endl;

        }

    out.close(); 

return 0;

delete[] RW;
delete[] RW_error;
}
#include <iostream>
#include <fstream>
#include <string>
#include "../../Librerie/random.h"
#include "../../Librerie/lib.h"

using namespace std;
 
int main (int argc, char *argv[]){

//Primo punto valutazione della media    

    int M=200000;
    int N= 100;
    int L= M/N;
    double* ave= new double[N]();
    double* ave2= new double[N]();
    double* err = new double [N]();
    double* sum_prog=new double [N]();

    Random rand;
    in_rand(rand);
 
        for (int i=0; i < N;i++){
            double appo=0;
            for(int j=0; j < L;j++){
                appo+=rand.Rannyu();
            }
            ave[i]=appo/L;
            ave2[i]=ave[i]*ave[i];
        }
        block_error(err,sum_prog,ave,ave2,N);

    ofstream out;
    out.open("Err.txt");
    
        for (int i=0; i<N;i++){
             out << err[i] << endl;
        }
    out.close();

    out.open("sum_prog.txt"); 

    for (int i=0; i<N;i++){
             out << sum_prog[i] << endl;
        }

    out.close();

    delete [] ave;
    delete [] ave2;
    delete [] err;
    delete [] sum_prog;
//--------------------------------------------------------------------------------------------------//

    //Secondo punto valutazione sigma^2
    double* var= new double[N]();
    double* var2= new double[N]();
    double* err_var = new double [N]();
    double* sum_var=new double [N]();

     for (int i=0; i < N;i++){
            double appo=0;
            
            for(int j=0; j < L;j++){
                appo += pow((rand.Rannyu()-0.5),2);
            }
            var[i]=appo/L;
            var2[i]=pow(var[i],2);
        }
    
     block_error(err_var,sum_var,var,var2,N);
    
    out.open("Err_var.txt");
    
        for (int i=0; i<N;i++){
             out << err_var[i] << endl;
        }
    out.close();

    out.open("sum_prog_var.txt"); 

    for (int i=0; i<N;i++){
             out << sum_var[i] << endl;
        }
    out.close();

    delete [] var;
    delete [] var2;
    delete [] err_var;
    delete [] sum_var;
//--------------------------------------------------------------------------------------------------//

    //Terzo punto chi quadro

    int M_c=100;
    int n_c=10000;
    int L_c = int(n_c/M_c);
    double* result=new double[L_c]();
    double appo_1=0;

     for (int i=0; i < L_c;i++){
        int* count= new int[M_c]();
        
        for(int j=0; j < n_c;j++){
                double k = rand.Rannyu()*100;
                int a = (int) k;
                count[a]=count[a]+1;
         }
        for(int k=0;k<M_c;k++){
            result[i] += pow(count[k]-L_c,2)/L_c;
             
        }
       // cout << result[i] << endl;
     }

    out.open("chi_q.txt");
    for(int k=0;k<M_c;k++){
        out << result[k] << endl;
       }

    delete [] result;

    out.close();


return 0;

}


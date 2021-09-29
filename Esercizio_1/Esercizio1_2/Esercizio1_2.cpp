#include <iostream>
#include <fstream>
#include <string>
#include "../../Librerie/random.h"
#include "../../Librerie/lib.h"

using namespace std;
 
int main (int argc, char *argv[]){

	int L = 10000;				//numero di realizzazioni
	int N1[4]={1,2,10,100};
	//inizio con N=1, ovvero un unico numero estratto secondo le diverse distribuzioni				//numero di punti che si considerano di volta in volta

	Random rand;
	in_rand(rand);
	//inizializzo il vettore degli sN relativi alla distribuzione uniforme a 0
	for(int h=0;h<4;h++){
		int N=N1[h];
		ofstream file_out;
		file_out.open("s"+to_string(N)+".txt");
		
		double *sN_u = new double[L]();
		double *sN_CL = new double[L]();
		double *sN_exp = new double[L]();
		
	
		for (int k=0; k<L; k++){
			
			//ciclo sugli N numeri casuali del file
			for (int j=0; j<N; j++){
				sN_u[k] += rand.Rannyu();
				sN_CL[k] += rand.cauchy(1,0);
				sN_exp[k] += rand.exp(1);
			}
			sN_u[k] /= N;
			sN_CL[k] /= N;
			sN_exp[k] /= N;
		}
		
		for (int k=0; k<L; k++){
			//esporto i valori delle L somme s2 per distribuzione uniforme, esponenziale e lorentziana in un file di testo "s2.txt"
			file_out << sN_u[k] << "    " << sN_exp[k] << "    " << sN_CL[k] << endl;
		}
		
		file_out.close();
		file_out.clear();




	delete[] sN_CL;
	delete[] sN_exp;
	delete[] sN_u;
}


    return 0;
}



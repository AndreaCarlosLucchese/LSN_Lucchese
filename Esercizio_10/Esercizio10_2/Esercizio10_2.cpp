#include <iostream>
#include <fstream>
#include <string>
#include "../../Librerie/random.h"
#include "../../Librerie/lib.h"
#include "TSP_SA.h"
#include <mpi.h>
using namespace std;
 


int main (int argc, char *argv[]){
    int nstep = 2000; //step totali dell'algoritmo
    int m_step= 1000; //step montecarlo
    int size, rank;
    int exchange, exchange2, exchange3;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status stat[4];
    MPI_Request req, req2;

    
   
    int itag=1;
    Random rand;
    in_rand(rand);

    SA_Algo square_tsp(32,1, "squared");
    ofstream out;
    out.open("outputs/city_square.dat");
    City cities_square = square_tsp.get_City();
    vector< array<double,2 > > locations_square=cities_square.get_Cities();
    for(int i=0; i<locations_square.size();i++)
        out << locations_square [i][0] << " " << locations_square[i][1] << endl;
    out.close();
    double fit1=square_tsp.getFitness();
    for(int istep = 0; istep < nstep; istep++){ 

	        square_tsp.Simluated_Annealing(m_step);
            
		    vector <unsigned int> bestpath = square_tsp.getPath().get_path(); //scrivo il best path 
            unsigned int* appo= &bestpath[0];
		    exchange = int(rand.Rannyu()*4);
		if(exchange == 1) {
            exchange2=2; 
            exchange3=3;
            }
		else if(exchange == 2) {
            exchange2=1; 
            exchange3=3;
            }
		else {
            exchange2=1; 
            exchange3=2;
            }

		if(rank == 0){ //faccio scambio tra rank 0 e rank==exchange
			MPI_Isend(appo, 32, MPI_INTEGER, exchange, itag, MPI_COMM_WORLD, &req);
			MPI_Recv(appo, 32, MPI_INTEGER, exchange, itag, MPI_COMM_WORLD, &stat[1]);
		        }
		else if(rank == exchange){
			MPI_Send(appo, 32, MPI_INTEGER, 0, itag, MPI_COMM_WORLD);
			MPI_Recv(appo, 32, MPI_INTEGER, 0, itag, MPI_COMM_WORLD, &stat[0]);
			}
		else if(rank == exchange2){ //faccio scambio tra rank==exchange2 e rank==exchange3
			MPI_Isend(appo, 32, MPI_INTEGER, exchange3, itag, MPI_COMM_WORLD, &req2);
			MPI_Recv(appo, 32, MPI_INTEGER, exchange3, itag, MPI_COMM_WORLD, &stat[3]);
			}
		else if(rank == exchange3){
			MPI_Send(appo, 32, MPI_INTEGER, exchange2, itag, MPI_COMM_WORLD);
			MPI_Recv(appo, 32, MPI_INTEGER, exchange2, itag, MPI_COMM_WORLD, &stat[2]);
			}

		square_tsp.ChangePath(bestpath); //sovrascrivo il bestpath scambiato
  }


   MPI_Finalize();
    out.open("outputs/square_new_best_path.dat");
    Salesman_Path new_best_square=square_tsp.getPath();
    vector<unsigned int> new_best_path_square = new_best_square.get_path();
    for(int i=0; i<new_best_path_square.size();i++)
        out << new_best_path_square[i] << endl;
    out.close();
    cout << "L1 migliore per quadrato post algoritmo: " << new_best_square.ComputeFit(cities_square)<< endl;
    cout << endl;

    out.open("outputs/square_step.dat");
    vector <double> L1 = square_tsp.get_L();
    for(int i=0; i<L1.size();i++){
        out << L1[i]<< endl;
    }
    
    return 0;
}

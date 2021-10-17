#include <iostream>
#include <fstream>
#include <string>
#include "../../Librerie/random.h"
#include "../../Librerie/lib.h"
#include "TSP_GA.h"
#include <mpi.h>
using namespace std;
 
// mpirun -np 4 ./Esercizio10_2.exe    comando da eseguire

int main (int argc, char *argv[]){
    int nstep = 100; //step totali dell'algoritmo
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

    Genetic_Algo square_tsp(400,32, "squared");
    ofstream out;
    out.open("outputs/city_square.dat");
    City cities_square = square_tsp.get_City();
    vector< array<double,2 > > locations_square=cities_square.get_Cities();
    for(int i=0; i<locations_square.size();i++)
        out << locations_square [i][0] << " " << locations_square[i][1] << endl;
    out.close();
    double fit1=square_tsp.getFitness();
    for(int istep = 0; istep < nstep; istep++){ 

	        square_tsp.Evolve();
            
		    vector <unsigned int> bestpath = square_tsp.getBestPath().get_path(); //scrivo il best path 
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

		square_tsp.changePath(bestpath); //sovrascrivo il bestpath scambiato
  }


   MPI_Finalize();
    out.open("outputs/square_new_best_path.dat");
    Salesman_Path new_best_square=square_tsp.getBestPath();
    vector<unsigned int> new_best_path_square = new_best_square.get_path();
    for(int i=0; i<new_best_path_square.size();i++)
        out << new_best_path_square[i] << endl;
    out.close();
    cout << "L1 migliore per quadrato post algoritmo: " << new_best_square.ComputeFit(cities_square)<< endl;
    cout << endl;

    
    return 0;
}

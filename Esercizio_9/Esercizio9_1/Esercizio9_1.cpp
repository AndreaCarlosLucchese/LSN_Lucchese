#include <iostream>
#include <fstream>
#include <string>
#include "../../Librerie/random.h"
#include "../../Librerie/lib.h"
#include "TSP_GA.h"
using namespace std;
 
void convergence(Genetic_Algo &, City &, string,string);


int main (int argc, char *argv[]){
    
    Genetic_Algo circle_tsp(400, 32, "circle");
    ofstream out;

    out.open("outputs/city_circle.dat");
    City cities = circle_tsp.get_City();
    vector< array<double,2 > > locations=cities.get_Cities();
    for(int i=0; i<locations.size();i++)
        out << locations[i][0] << " " << locations[i][1] << endl;
    out.close();

    Salesman_Path best=circle_tsp.getBestPath();

    out.open("outputs/circle_initial_best_path.dat");
    vector<unsigned int> best_path = best.get_path();
    for(int i=0; i<best_path.size();i++)
        out << best_path[i] << endl;
    out.close();
    cout << "L1 migliore per circonferenza all'inizio: " << best.ComputeFit(cities)<< endl;
    cout << endl;

    convergence(circle_tsp,cities,"outputs/av_circle.dat","outputs/best_circle.dat");
    out.open("outputs/circle_new_best_path.dat");
    Salesman_Path new_best=circle_tsp.getBestPath();
    vector<unsigned int> new_best_path = new_best.get_path();
    for(int i=0; i<new_best_path.size();i++)
        out << new_best_path[i] << endl;
    out.close();
    cout << "L1 migliore per circonferenza post algoritmo: " << new_best.ComputeFit(cities)<< endl;
    cout << endl;


    //Square 32 cities 

    Genetic_Algo square_tsp(400, 32, "squared");
    out.open("outputs/city_square.dat");
    City cities_square = square_tsp.get_City();
    vector< array<double,2 > > locations_square=cities_square.get_Cities();
    for(int i=0; i<locations_square.size();i++)
        out << locations_square [i][0] << " " << locations_square[i][1] << endl;
    out.close();

    Salesman_Path best_square=square_tsp.getBestPath();

    out.open("outputs/square_initial_best_path.dat");
    vector<unsigned int> best_path_square = best_square.get_path();
    for(int i=0; i<best_path_square.size();i++)
        out << best_path_square[i] << endl;
    out.close();
    cout << "L1 migliore per quadrato all'inizio: " << best_square.ComputeFit(cities_square)<< endl;
    cout << endl;

    convergence(square_tsp,cities,"outputs/av_square.dat","outputs/best_square.dat");
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

void convergence(Genetic_Algo& tsp, City& cities,string out,string out1 ){ // algoritmo di convergenza 
    double fit1, fit2;
    ofstream av,best;
    av.open(out);
    best.open(out1);
    fit1=tsp.getBestPath().ComputeFit(cities);
    for(unsigned int i = 0; i<=5000; i++){
        tsp.Evolve();
        if(i%100==0 and i!=0){
            fit2=tsp.getBestPath().ComputeFit(cities);
            cout << "Numero di iterazioni :" << i << endl;
            if(fit1-fit2< 0.0000001){
                cout << "Convergenza dopo : " << i <<" iterazioni"<<endl;
                break;
            }
            fit1=fit2;
        }
        av<<tsp.getAvFitness()<< endl;
        best <<tsp.getFitness()<<endl;
        if(i==5000){
            cout << "L'algoritmo non coverge " << i <<endl;
            exit(0);
        }
        
    }

    av.close();
    best.close();
}

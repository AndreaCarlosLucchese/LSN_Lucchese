#include <iostream>
#include <fstream>
#include <string>
#include "../../Librerie/random.h"
#include "../../Librerie/lib.h"
#include "TSP_SA.h"
using namespace std;
 
//void convergence(Genetic_Algo &, City &, string);


int main (int argc, char *argv[]){
    
    SA_Algo circle_tsp(32,1, "circle");
    ofstream out;

    out.open("outputs/city_circle.dat");
    City cities = circle_tsp.get_City();
    vector< array<double,2 > > locations=cities.get_Cities();
    for(int i=0; i<locations.size();i++)
        out << locations[i][0] << " " << locations[i][1] << endl;
    out.close();

    Salesman_Path best=circle_tsp.getPath();

    out.open("outputs/circle_initial_best_path.dat");
    vector<unsigned int> best_path = best.get_path();
    for(int i=0; i<best_path.size();i++)
        out << best_path[i] << endl;
    out.close();
    cout << "L1 migliore per circonferenza all'inizio: " << best.ComputeFit(cities)<< endl;
    cout << endl;

    
    circle_tsp.Simluated_Annealing();

    out.open("outputs/circle_new_best_path.dat");
    Salesman_Path new_best=circle_tsp.getPath();
    vector<unsigned int> new_best_path = new_best.get_path();
    for(int i=0; i<new_best_path.size();i++)
        out << new_best_path[i] << endl;
    out.close();
    cout << "L1 migliore per circonferenza post algoritmo: " << new_best.ComputeFit(cities)<< endl;
    cout << endl;
    out.open("outputs/circle_step.dat");
    vector <double> L1 = circle_tsp.get_L();
    for(int i=0; i<L1.size();i++){
        out << L1[i]<< endl;
    }
    L1.clear();
    out.close();

    //Square 32 cities 

    SA_Algo square_tsp(32,1, "squared");
    out.open("outputs/city_square.dat");
    City cities_square = square_tsp.get_City();
    vector< array<double,2 > > locations_square=cities_square.get_Cities();
    for(int i=0; i<locations_square.size();i++)
        out << locations_square [i][0] << " " << locations_square[i][1] << endl;
    out.close();

    Salesman_Path best_square=square_tsp.getPath();
    out.open("outputs/square_initial_best_path.dat");
    vector<unsigned int> best_path_square = best_square.get_path();
    for(int i=0; i<best_path_square.size();i++)
        out << best_path_square[i] << endl;
    out.close();
    cout << "L1 migliore per quadrato all'inizio: " << best_square.ComputeFit(cities_square)<< endl;
    cout << endl;

    square_tsp.Simluated_Annealing();
    out.open("outputs/square_new_best_path.dat");
    Salesman_Path new_best_square=square_tsp.getPath();
    vector<unsigned int> new_best_path_square = new_best_square.get_path();
    for(int i=0; i<new_best_path_square.size();i++)
        out << new_best_path_square[i] << endl;
    out.close();
    cout << "L1 migliore per quadrato post algoritmo: " << new_best_square.ComputeFit(cities_square)<< endl;
    cout << endl;

    out.open("outputs/square_step.dat");
    L1 = square_tsp.get_L();
    for(int i=0; i<L1.size();i++){
        out << L1[i]<< endl;
    }
    

    return 0;
}
#include <cmath>
#include <cstdlib>
#include <string>
#include "../../Librerie/random.h"
#include "../../Librerie/lib.h"
#include <algorithm>
#include <vector>
#include <array>
#include <iterator>
#ifndef _TSP_GA_H_
#define _TSP_GA_H_
using namespace std;

class City{
    private:
        vector< array<double,2 > > m_city;
        Random m_rand;

    public:
        City(unsigned int, string);
        double get_L1(unsigned int, unsigned int);
        vector<array <double,2> > get_Cities();

};


class Salesman_Path{

    private:
        vector<unsigned int> m_path;
    
        double m_fitness;
        Random m_rand;
    public:
        Salesman_Path(unsigned int);
        Salesman_Path(vector<unsigned int>);
        double ComputeFit(City);
        vector<unsigned int>  get_path();
        void Mutation(Random rand);
        void Single_permutation(Random rand);
        void Shift(Random rand);
        void Inversion(Random rand);
        bool Check()const;
};


class Genetic_Algo{
    private:
        vector <Salesman_Path> m_population;
        vector <Salesman_Path> new_population;
        Random m_rand;
        City cities;
    public:
        Genetic_Algo(unsigned int , unsigned int ,string);
        unsigned int Selection();
        void Evolve();
        City get_City();
        Salesman_Path getBestPath();
        double getFitness();
        double getAvFitness();
};


#endif

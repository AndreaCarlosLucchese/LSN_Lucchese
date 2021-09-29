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
        void Mutation(Random );
        void Single_permutation(Random );
        void Shift(Random );
        void Inversion(Random );
        bool Check()const;
};

class SA_Algo{
    private:
        Salesman_Path path;
        vector <double> L1;
        double m_T;
        Random m_rand;
        City cities;
    public:
        SA_Algo(unsigned int,double ,string);
        Salesman_Path getPath();
        double getFitness();
        void Simluated_Annealing(unsigned int);
        City get_City();
        vector <double> get_L();
        void ChangePath(vector <unsigned int>);
};



#endif

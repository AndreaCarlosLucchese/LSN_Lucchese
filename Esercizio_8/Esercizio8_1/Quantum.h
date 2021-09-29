#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include "../../Librerie/random.h"
#include "../../Librerie/lib.h"
#include <algorithm>
#ifndef __Quantum_hole__
#define __Quantum_hole__

class Quantum_hole {

    private:

    double m_x;
    double m_mu, m_sigma,m_step;
    Random m_rand;
    int m_counter;

    public:
    Quantum_hole(double x, double mu, double sigma, double step);
    ~Quantum_hole();
    double get_radius();
    void reset_start(double x,double mu, double sigma, double step);
    double prob(double x);
    double en(double x);
    void metropolis_unif();
    int get_counter();
    double get_step();
};

#endif 
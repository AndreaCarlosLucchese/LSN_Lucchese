#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include "../../Librerie/random.h"
#include "../../Librerie/lib.h"
#include <algorithm>
#ifndef __hydro__
#define __hydro__

class Hydro {

    private:

    double m_x,m_y,m_z, m_step,m_p;
    Random m_rand;
    int m_counter;

    public:
    Hydro(double x, double y, double z, double step, double p);
    ~Hydro();
    double get_radius();
    void reset_start(double x, double y, double z,double step);
    void metropolis_unif();
    void metropolis_gauss();
    double prob_100(double x, double y, double z);
    double prob_210(double x, double y, double z);
    int get_counter();
    double get_step();
};

#endif 

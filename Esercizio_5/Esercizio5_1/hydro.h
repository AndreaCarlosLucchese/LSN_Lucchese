#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <vector>
#include "../../Librerie/random.h"
#include "../../Librerie/lib.h"
#include <algorithm>
#include <iostream>
#ifndef __hydro__
#define __hydro__
using namespace std;
class Hydro {

    private:
    vector <double> m_x,m_y,m_z;
    vector <double> m_r;
    double m_step,m_p;
    Random m_rand;
    int m_counter;

    public:
    Hydro(double x, double y, double z, double step, double p);
    ~Hydro();
    double get_radius(int);
    void reset_start(double x, double y, double z,double step);
    vector <double>  Equilibration(int);
    void metropolis_unif(int);
    void metropolis_gauss(int);
    double prob_100(double x, double y, double z);
    double prob_210(double x, double y, double z);
    int get_counter();
    double get_step();
    double* get_pos();
};

#endif 

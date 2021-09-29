#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include "random.h"

//dist

double dist(double x1, double x2, double y1, double y2);

// block error

void block_error(double* error, double* sum_prog,double* x, double* x2, int n);


//Inzializzatore random

void in_rand(Random &rnd);

// random walk discreto

void RWD(double* RW, double* RW_error, Random &rnd, int f);

// random walk continuo

void RWC(double* RW, double* RW_error, Random &rnd, int f);

// d1

double d1(double S, double T, double t,double K, double r, double sigma);

// d2
double d2(double d1, double T, double t, double sigma);

// N(x)
double N(double x);

// Call-option
double C(double S, double T, double t, double r, double N1, double N2, double K);

//put-option

double P(double S, double T, double t, double r, double N1, double N2, double K);


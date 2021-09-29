#include "Quantum.h"

using namespace std;

Quantum_hole::Quantum_hole(double x, double mu, double sigma, double step){

        m_x=x;
        m_mu=mu;
        m_sigma=sigma;
        m_step = step;
        in_rand(m_rand);
        m_counter=0;
}

Quantum_hole::~Quantum_hole(){

}

double Quantum_hole::get_radius(){

    return m_x;
}

void Quantum_hole::reset_start(double x, double mu, double sigma,double step){
    m_x = x;
    m_step=step;
    m_mu=mu;
    m_sigma=sigma;
    m_counter=0;
   // in_rand(m_rand);
}

double Quantum_hole::prob(double x){ // distribuzione di probabilità della funzione d'onda
   double exp_1=exp(-pow(x+m_mu,2)/(2*m_sigma*m_sigma))+exp(-pow(x-m_mu,2)/(2*m_sigma*m_sigma));
    return pow(exp_1,2);

}

void Quantum_hole::metropolis_unif(){ //metropolis
    double x_new,a;
        x_new=m_x+(m_rand.Rannyu()-0.5)*m_step;
        a=min(1.,prob(x_new)/prob(m_x));
        if(m_rand.Rannyu()<a){
            m_x=x_new;
            m_counter++;
        
    }

}

double Quantum_hole::en(double x){ // calcolo dell'eneriga in funzione di una posizione 
    double exp1=exp(-pow(x+m_mu,2)/(2*m_sigma*m_sigma));
    double exp2=exp(-pow(x-m_mu,2)/(2*m_sigma*m_sigma));
    double appo= (1./(2. * pow(m_sigma,2))) * (1. - (pow(x,2) + pow(m_mu,2))/pow(m_sigma,2) - (2 * x * m_mu/pow(m_sigma,2)) * (exp1 - exp2)/(exp1+ exp2));
    return appo+pow(x,4)-pow(x,2)*2.5;
}

int Quantum_hole::get_counter(){
    return m_counter;
}

double Quantum_hole::get_step(){  
    return m_step;
}




#include "hydro.h"

using namespace std;

Hydro::Hydro(double x, double y, double z, double step, double p){

        m_x=x;
        m_y=y;
        m_z=z;
        m_step = step;
        in_rand(m_rand);
        m_p=p;
        m_counter=0;
}

Hydro::~Hydro(){

}

double Hydro::get_radius(){

    return sqrt(m_x*m_x + m_y*m_y + m_z*m_z);
}

void Hydro::reset_start(double x, double y, double z,double step){
    m_x = x;
    m_y = y;
    m_z = z;
    m_step=step;
    m_counter=0;
   // in_rand(m_rand);
}

double Hydro::prob_100(double x, double y, double z){ //funzione d'onda ground state
    double radius=sqrt(x*x + y*y+ z*z);
    return exp(-2.*radius)/M_PI;
}

double Hydro::prob_210(double x, double y, double z){ //funzione d'onda primo stato eccitato
    double radius=sqrt(x*x + y*y + z*z);
    return exp(-radius)*x*x/(32.*M_PI);
}

void Hydro::metropolis_unif(){ //algoritmo metropolis con distribuzione uniforme
    double x_new,y_new,z_new,a;
    if(m_p==0){
        x_new=m_x+(m_rand.Rannyu()-0.5)*m_step;
        y_new=m_y+(m_rand.Rannyu()-0.5)*m_step;
        z_new=m_z+(m_rand.Rannyu()-0.5)*m_step;

        a=min(1.,prob_100(x_new,y_new,z_new)/prob_100(m_x,m_y,m_z));
        if(m_rand.Rannyu()<a){
            m_x=x_new;
            m_y=y_new;
            m_z=z_new;
            m_counter++;
        }
    }
    else{
        x_new=m_x+(m_rand.Rannyu()-0.5)*m_step;
        y_new=m_y+(m_rand.Rannyu()-0.5)*m_step;
        z_new=m_z+(m_rand.Rannyu()-0.5)*m_step;

        a=min(1.,prob_210(x_new,y_new,z_new)/prob_210(m_x,m_y,m_z));
        if(m_rand.Rannyu()<a){
            m_x=x_new;
            m_y=y_new;
            m_z=z_new;
            m_counter++;
        }

    }

}


void Hydro::metropolis_gauss(){ //algoritmo metropolis con distribuzione gaussiana
    double x_new,y_new,z_new,a;
    if(m_p==0){
        x_new=m_x+m_rand.Gauss(0,1)*m_step;
        y_new=m_y+m_rand.Gauss(0,1)*m_step;
        z_new=m_z+m_rand.Gauss(0,1)*m_step;

        a=min(1.,prob_100(x_new,y_new,z_new)/prob_100(m_x,m_y,m_z));
        if(m_rand.Rannyu()<a){
            m_x=x_new;
            m_y=y_new;
            m_z=z_new;
            m_counter++;
        }

    }
    else{

        x_new=m_x+m_rand.Gauss(0,1)*m_step;
        y_new=m_y+m_rand.Gauss(0,1)*m_step;
        z_new=m_z+m_rand.Gauss(0,1)*m_step;

        a=std::min(1.,prob_210(x_new,y_new,z_new)/prob_210(m_x,m_y,m_z));
        if(m_rand.Rannyu()<a){
            m_x=x_new;
            m_y=y_new;
            m_z=z_new;
            m_counter++;
        }

    }
}


int Hydro::get_counter(){
    return m_counter;
}

double Hydro::get_step(){  
    return m_step;
}











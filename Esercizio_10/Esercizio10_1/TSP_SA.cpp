#include "TSP_SA.h"
using namespace std;


//City location

City::City(unsigned int n_cities, string location): m_city(n_cities){

    in_rand(m_rand);    

    if(location == "circle"){//raggio=1
        double theta;
        for (unsigned int i =0; i<n_cities; ++i){
            theta=m_rand.Rannyu()*M_PI*2.;
            m_city[i][0]=cos(theta);
            m_city[i][1]=sin(theta);
        }
    }
    else if (location == "squared"){
        for (unsigned int i =0; i<n_cities; ++i){
            m_city[i][0]=m_rand.Rannyu();
            m_city[i][1]=m_rand.Rannyu();
        }
    }



}

vector< array <double , 2 > > City::get_Cities() {
    return m_city;
}

double City::get_L1(unsigned int a, unsigned int b){
    double x = m_city[a][0]-m_city[b][0];
    double y = m_city[a][1]-m_city[b][1];
    return sqrt(x*x+y*y);
}

//SalesmanPath

Salesman_Path::Salesman_Path(unsigned int n_cities):m_path(n_cities){
    in_rand(m_rand);
    for(int i=0; i<n_cities;i++){
        m_path[i]=i;
    }
    random_shuffle(m_path.begin() ,m_path.end());
    
    m_fitness=0.;
}

Salesman_Path::Salesman_Path(vector<unsigned int> path){
    m_path=path;
}

double Salesman_Path::ComputeFit(City cities){
        m_fitness=0.;
        for(unsigned int i=0; i<m_path.size()-1; ++i)
            m_fitness+=cities.get_L1(m_path[i], m_path[i+1]);
        m_fitness+=cities.get_L1(m_path[m_path.size()-1], m_path[0]); 
    return m_fitness;
}

vector <unsigned int> Salesman_Path::get_path(){
    return m_path;
}

void Salesman_Path::Mutation(Random rand){
    unsigned int index=int(rand.Rannyu()*3);
    if (index==0)
            Single_permutation(rand);
    else if(index==1)
            Shift(rand);
    else if(index==2)
            Inversion(rand);
}

void Salesman_Path::Single_permutation(Random rand){
    unsigned int temp;
    unsigned int i=int(rand.Rannyu()*m_path.size());
    unsigned int j=int(rand.Rannyu()*m_path.size());
    do
    {
        j=int(rand.Rannyu()*m_path.size());
    } while (i!=j);
    temp=m_path[i];
    m_path[i]=m_path[j];
    m_path[j]=temp;
}

void Salesman_Path::Shift(Random rand){
    unsigned int pos=int(rand.Rannyu()*(m_path.size()-1))+1; 
    rotate ( m_path.begin(), m_path.begin()+pos, m_path.end() ); 
}

void Salesman_Path::Inversion(Random rand){
    unsigned int n=int(rand.Rannyu()*(m_path.size()-1))+2; 
    unsigned int pos=int(rand.Rannyu()*(m_path.size()-n));
    reverse(m_path.begin()+pos, m_path.begin()+pos+n);
}


bool Salesman_Path::Check() const {
    unsigned int tot=0;
    unsigned int check=0;
    for(int i=0; i<m_path.size();i++){
        unsigned int count=0;
        tot+=m_path[i];
        check++;
        for(int j=0; j<m_path.size();j++){
            if(m_path[i]==j)
                count++;
            if(count >1)
                return false;
        }
    }
    if(tot!=check)
        return false;
    return true;
}

//Simulated Annealing Algorithm

SA_Algo::SA_Algo( unsigned int ncities, double T, string loc): path(ncities), cities(ncities,loc) {
        in_rand(m_rand);
        m_T=T;

}


void SA_Algo::Simluated_Annealing(){ // Algoritmo di Simulated Annealing con annessa convergenza
    double prob=0;
    for(int i=1; i<=100000;i++){

        double appo = double(m_T/pow(10,i));
        cout << appo << endl;
        double fit1=getFitness();

        for(int j=0; j< 1000; j++){
        Salesman_Path old_path=path;
        path.Mutation(m_rand);
        prob = exp(-(1/appo)*(path.ComputeFit(cities)-old_path.ComputeFit(cities)));
        if(prob < m_rand.Rannyu())
            path=old_path;
        L1.push_back(path.ComputeFit(cities));
        }
        double fit2=getFitness();
        if(abs(fit2-fit1) < 0.000001){
            cout << "Convergenza raggiunta alla temperatura di " << appo << endl;
            cout << endl;
            break;
        }
    }
}


City SA_Algo::get_City(){
    return cities;
}

Salesman_Path SA_Algo::getPath(){
   return path;
}

double SA_Algo::getFitness(){
    return path.ComputeFit(cities);
}

vector <double> SA_Algo::get_L(){
    return L1;
}

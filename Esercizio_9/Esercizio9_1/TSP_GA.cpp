#include "TSP_GA.h"
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

vector< array <double , 2 > > City::get_Cities() { // mi da la posizione di tutte le città
    return m_city;
}

double City::get_L1(unsigned int a, unsigned int b){ // calcola la funzione costo
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
    random_shuffle(m_path.begin()+1 ,m_path.end());
    
    m_fitness=0.;
}

Salesman_Path::Salesman_Path(vector<unsigned int> path){ 
    m_path=path;
}

double Salesman_Path::ComputeFit(City cities){ // calcolo il fit del path
        m_fitness=0.;
        for(unsigned int i=0; i<m_path.size()-1; ++i)
            m_fitness+=cities.get_L1(m_path[i], m_path[i+1]);
        m_fitness+=cities.get_L1(m_path[m_path.size()-1], m_path[0]); 
    return m_fitness;
}

vector <unsigned int> Salesman_Path::get_path(){ // restituisce il path di un singolo salesman
    return m_path;
}

void Salesman_Path::Mutation(Random rand){ // Mutazioni possibili
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


bool Salesman_Path::Check() const { // Funzione per checkare la validità di un path
    unsigned int tot=0;
    unsigned int check=0;
    for(int i=0; i<m_path.size();i++){
        unsigned int count=0;
        tot+=m_path[i];
        check+=i;
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

//Genetic Algorithm

Genetic_Algo::Genetic_Algo(unsigned int npopulation,unsigned int ncities, string loc): cities(ncities,loc){
        in_rand(m_rand);
        for(int i=0; i<npopulation;i++){
            m_population.push_back(Salesman_Path(ncities));
            if(m_population[i].Check()==false){
            cout << "Errore nella generazione" << endl;
            exit(0);
            }
        }
        for(int i=0;i<npopulation;i++){  //ordino le città in funzione del costo
            for(int j=i+1;j<npopulation;j++){
                if(m_population[i].ComputeFit(cities) > m_population[j].ComputeFit(cities))
                    swap(m_population[i],m_population[j]);
            }
        }
    }

unsigned int Genetic_Algo::Selection(){ // Selezione 
        return m_population.size()*pow(m_rand.Rannyu(),5);
        
}

void Genetic_Algo::Evolve(){

    for(int i=0; i< m_population.size()/2;i++){  //Crossover
        unsigned int sel1=Selection();
        unsigned int sel2=Selection();
        if(m_rand.Rannyu()<0.5){
            vector <unsigned int> p1=m_population[sel1].get_path();
            vector <unsigned int> p2=m_population[sel2].get_path();
            unsigned int index=int(m_rand.Rannyu()*p1.size());
            vector <unsigned int > index1;
            vector <unsigned int > index2;
            vector<unsigned int>::iterator it;
            for(int i=index; i<p1.size();i++){
                it=find(p2.begin(),p2.end(),p1[i]);
                index1.push_back(it-p2.begin());
                it=find(p1.begin(),p1.end(),p2[i]);
                index2.push_back(it-p1.begin());
            }

            sort(index1.begin(),index1.end());
            sort(index2.begin(),index2.end());
            vector <unsigned int> new_p1,new_p2;
            new_p1=p1;
            new_p2=p2;
            for(int i=index;i<p1.size();i++){
                new_p1[i]=p2[index1[i-index]];
                new_p2[i]=p1[index2[i-index]];
            }

            new_population.push_back(Salesman_Path(new_p1));
            new_population.push_back(Salesman_Path(new_p2));
            if(new_population[i].Check()==false){
                cout << "Errore nel crossover" << endl;
                exit(0);
            }
            if(new_population[i+1].Check()==false){
                cout << "Errore nel crossover" << endl;
                exit(0);
            }
        }
        else{
            new_population.push_back(m_population[sel1]);
            new_population.push_back(m_population[sel2]);
        }
    }

    for(int i=0; i< new_population.size();i++){ //Mutazione
        if(m_rand.Rannyu()<0.1){
            new_population[i].Mutation(m_rand);
        }
    }

    for(int i=0;i<new_population.size();i++){    //Rivalutazione di L
            for(int j=i+1;j<new_population.size();j++){
                if(new_population[i].ComputeFit(cities) > new_population[j].ComputeFit(cities))
                    swap(new_population[i],new_population[j]);
            }
        }   
   m_population=new_population;
   new_population.clear();
}

Salesman_Path Genetic_Algo::getBestPath(){
    return m_population[0];
}

double Genetic_Algo::getFitness(){
    return m_population[0].ComputeFit(cities);
}

City Genetic_Algo::get_City(){
    return cities;
}

double Genetic_Algo::getAvFitness(){ // Fitness di tutta la popolazione
    double sum=0;
    for(int i=0; i<m_population.size()/2;i++){
        sum+=m_population[i].ComputeFit(cities);
    }
    return sum/(m_population.size()/2);
}
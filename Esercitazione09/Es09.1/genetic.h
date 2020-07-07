#ifndef __Genetic__
#define __Genetic__


#include <vector>
#include "random.h"


using namespace std;


void InitRandom(Random *rnd);

class Individuo {

private:

	vector<int> path;
	int N;
	double fitness;

	
protected:

public:

	// constructors
	Individuo();				//vettore [1] (fitness=0)
	Individuo(int N);			//vettore [1, 2, ..., N] (fitness=0)
	Individuo(vector<int> v);		//vettore passato in input (fitness=0)
	Individuo(const Individuo &ind);	//copy constructor

	// destructor
	~Individuo();

	Individuo& operator=(const Individuo&);
	bool operator<(const  Individuo & ind) { return fitness < ind.get_fitness(); }

	// methods



	bool check();				//true se rispetta le condizioni, false altrimenti
	void reset();				//resetta al vettore [1, 2, ..., N]
	void reset(vector<int> v);		//resetta al vettore v passato
	void reset(Individuo ind);		//resetta all'individuo ind passato

	void printInd() const;

	vector<int> get_vec() const{ return path; }
	int get_comp(int i) const{ return path[i]; }
	int get_size() const{ return N; };
	double get_fitness() const{ return fitness; };
	void set_fitness(vector<double> xc, vector<double> yc){ fitness = L1(xc, yc);}		//scegliere se usare L1 o L2

	void perm(int i, int j);			//permutazione dell'elemento in posizione i e quello in posizione j
	void shift(int n, int i, int j);		//shift degli elementi che vanno da i a j di n posizioni
	void perm_m(int m, int i, int j);		//permutazione degli m elementi dalla posizione i con gli m elementi dalla posizione j
	void inversion(int i, int j);			//inversione degli elementi compresi tra i e j

	bool mutate(double prob, Random *rnd);		//prova ciascuno dei 4 operatori con probabilità prob

	double L1(vector<double> xc, vector<double> yc) const;
	double L2(vector<double> xc, vector<double> yc) const;	

};



class Popolazione{

private:
 
	int I;		//numero di individui
	int N;		//lunghezza di ciascun individuo

	vector<Individuo> pop;
	vector<Individuo> new_gen;
	vector<double> xc;		//coordinate x delle città
	vector<double> yc;		//coordinate y delle città

	bool sorted;	//true se il vettore è ordinato
	bool updated;	//true se tutte le fitness sono aggiornate

	Random rnd;

protected:

public:
	// constructors
	Popolazione();				//parte da [1, 2, ..., N] e crea la popolazioni con permutazioni casuali successive (I=100, N=32)
	Popolazione(int Is, int Ns);		// ''			''			''			''	    (I e N passato in input)
	Popolazione(Individuo ind);		//parte dall'individuo passato e crea la popolazioni con permutazioni casuali successive (I=100, N=32)
	Popolazione(int Is, Individuo ind);	// ''			''			''			''	    	 (I passato in input)
	Popolazione(const Popolazione &p);	//copy constructor

	// destructor
	~Popolazione();


	Popolazione& operator=(const Popolazione& p);

	// methods

	void set_cities(vector<double> x, vector<double> y);
	void set_random(Random r){ rnd=r; };

	int get_I() const{ return I; };						//restituisce il numero di individui
	Individuo get_ind(int i) const;						//restituisce l'i-esimo indivudo (solo lettura)
	double get_fitness(int i) const{ return pop[i].get_fitness();}		//fitness dell'i-esimo individuo
	vector<double> get_xc() const{ return xc; };				//coordinate x delle città
	vector<double> get_yc() const{ return yc; };				//coordinate y delle città


	bool check_ind(int i);			//true se l'i-esimo individuo rispetta le condizione, false altrimenti
	bool IsSorted()const{ return sorted;};
	void printPop()const;
	void printNG()const;

	void update_fitness();
	void update_fitness(int i);
	void update_fitness_ng();
	void update_fitness_ng(int i);

	void sort_pop();
	void quicksort(int primo, int ultimo);	//ordina gli individuo per fitness, dalla più bassa (migliore) alla più alta (peggiore)
	void scambia(int i, int j);		//scambia individui i e j
	int selection_unif();
	int selection();			//seleziona un indice di un individuo in base alla fitness


	void mutate_ind(int i, double prob_m);			//tenta di mutare l'i-esimo individuo, con probabilità prob_m
	void mutate_all(double prob_m);				//tenta di mutare tutti gli individui
	void crossover(int i, int j, double prob_c);		//crossover tra l'individuo i-esimo e l'individuo j-esimo, con probabilità prob_c
	void update_gen();					//sostituisce la vecchia generazione con la nuova
	
	void RS(double prob_m);			//random search
	void GA(double prob_m, double prob_c);	//algoritmo genetico


};

#endif // __Genetic__

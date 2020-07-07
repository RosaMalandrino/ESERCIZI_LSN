#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
//#include <vector>
#include <algorithm>
#include "genetic.h"

using namespace std;



void InitRandom(Random *rnd){

	int seed[4];
	int p1, p2;
	ifstream Primes("Primes");
	if (Primes.is_open()){
		Primes >> p1 >> p2 ;
	} else cerr << "PROBLEM: Unable to open Primes" << endl;
	Primes.close();

	ifstream input("seed.in");
	string property;
	if (input.is_open()){
		while ( !input.eof() ){
         		input >> property;
         		if( property == "RANDOMSEED" ){
            			input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            			rnd->SetRandom(seed,p1,p2);
         		}
      		}
      	input.close();
   	} else cerr << "PROBLEM: Unable to open seed.in" << endl;

}

/***************************************** COSTRUTTORI E DISTRUTTORE ********************************************/

Individuo :: Individuo(){

	N=1;
	vector<int> v(N,1);

	if(!path.empty()){
		path.clear();
		path.resize(0);
		path.shrink_to_fit();
	}
	path=v;

	fitness=0.;


}

Individuo :: Individuo(int dim){

	N=dim;
	vector<int> v(N);
	for(int i=0; i<N; i++)	v[i]=i+1;

	if(!path.empty()){
		path.clear();
		path.resize(0);
		path.shrink_to_fit();
	}
	path=v;

	fitness=0.;


}

Individuo :: Individuo(vector<int> v){

	N=v.size();

	if(!path.empty()){
		path.clear();
		path.resize(0);
		path.shrink_to_fit();
	}
	path=v;

	fitness=0.;


}

Individuo :: Individuo(const Individuo &ind){

	N=ind.get_size();

	if(!path.empty()){
		path.clear();
		path.resize(0);
		path.shrink_to_fit();
	}
	path=ind.get_vec();

	fitness=ind.fitness;

}


Individuo :: ~Individuo(){}


Individuo& Individuo::operator=(const Individuo& ind){

	N = ind.N;
	
	if(!path.empty()){
		path.clear();
		path.resize(0);
		path.shrink_to_fit();
	}
	path=ind.path;

	fitness = ind.fitness;

	
	return *this;

}




/**************************************************** METODI ***************************************************************/

bool Individuo :: check(){

	bool ok=true;

	for(int i=0; i<N; i++){


		if(path[i]<1 or path[i]>N){ ok=false; break;}

		for(int j=i+1; j<N; j++){

			if(path[j]==path[i]){ ok=false; break;}

		}

	}


	if(ok==true and path[0]!=1){

		//cout<<"Rotazione del percorso finché 1 non è il primo elemento"<<endl;

		while(path[0]!=1){		//se il vettore va bene, ma 1 non è il primo elemento

			int a=path[0];
			path.push_back(a);
			path.erase(path.begin());
			
		}
	
	}

	return ok;

}



void Individuo :: reset(){

	vector<int> v(N);
	for(int i=0; i<N; i++)	v[i]=i+1;
	path=v;	

}

void Individuo :: reset(vector<int> v){

	N=v.size();
	path=v;	

}


void Individuo :: reset(Individuo ind){

	N=ind.get_size();
	path=ind.get_vec();	

}


void Individuo :: printInd() const{

	cout<<"Path: [ ";
	for(int i=0; i<N-1; i++)	cout<<path[i]<<", ";
	cout<<path[N-1]<<" ] Fitness: "<<fitness<<endl;
}


void Individuo :: perm(int i, int j){

	if(i<=0 or i>=N or j<=0 or j>=N){
		cerr<<"Gli indici "<<i<<", "<<j<<" sono fuori dal range consentito, permutazione non eseguita"<<endl;
		return;
	}

	int a=path[i];
	path[i]=path[j];
	path[j]=a;


}


void Individuo :: shift(int n, int i, int j){

	if(n<0){
		cerr<<"Shift ammessi n>=0"<<endl;
		return;
	}

	if(i<=0 or i>=N or j<=0 or j>=N){
		cerr<<"Gli indici "<<i<<", "<<j<<" sono fuori dal range consentito, permutazione non eseguita"<<endl;
		return;
	}

	int low=i, high=j;
	if(j<i){ low=j; high=i; }	//controllo ordine

	//int m = high - low +1;		//numero di elementi da shiftare

	if(high+n>=N){
		cerr<<"L'indice superiore "<<high<<" non è seguito da abbastanza elementi per shiftare un blocco di "<<n<<" elementi: ";
		cerr<<"shift non eseguito"<<endl;
		return;
	}

	
	for(int k=0; k<n; k++){

		int a = path[high+k+1];
		path.erase(path.begin() + high+k+1);
		path.insert(path.begin() + low + k, a);
	
	}

}



void Individuo :: perm_m(int m, int i, int j){

	if(m>(N-1)/2){
		cerr<<"Blocco troppo grande per le dimensioni del vettore, impossibile permutare"<<endl;
		return;
	}

	if(i<=0 or i>=N or j<=0 or j>=N){
		cerr<<"Gli indici "<<i<<", "<<j<<" sono fuori dal range consentito, permutazione non eseguita"<<endl;
		return;
	}

	int low=i, high=j;
	if(j<i){ low=j; high=i; }	//controllo ordine

	if(high+m>N){
		cerr<<"L'indice superiore "<<high<<" non è seguito da abbastanza elementi per permutare un blocco di "<<m<<" elementi: ";
		cerr<<"permutazione non eseguita"<<endl;
		return;
	}

	if(high-low<m){
		cerr<<"Gli indici "<<low<<", "<<high<<" sono troppo vicini per permutare un blocco di "<<m<<" elementi: ";
		cerr<<"permutazione non eseguita"<<endl;
		return;
	}

	for(int k=0; k<m; k++){

		perm(low+k, high+k);
	
	}

}


void Individuo :: inversion(int i, int j){

	int low=i, high=j;	
	if(j<i){ low=j; high=i; }	//controllo ordine 

	if(low<=0)	low=1;		//se low è fuori dal range, tenendo conto che il primo elemento deve rimanere invariato
	if(high>=N)	high=N-1;	//se high è fuori dal range pongo high all'ultimo elemento
					//se gli indici sono fuori range viene invertito tutto il vettore eccetto la posizione 0

	int m = high - low +1;		//numero di elementi da permutare

	for(int k=0; k<m/2; k++){

		perm(low+k, high-k);

	}


}



bool Individuo :: mutate(double prob_m, Random *rnd){


	vector<int> old(N);
	old=get_vec();
	bool mutated = false;

	if(rnd->Rannyu()<prob_m){
		
		int j = 1 + int(rnd->Rannyu(0, N-1));	//intero tra 1 e N-1
		int k;

		do{ k = 1 + int(rnd->Rannyu(0, N-1));	//intero tra 1 e N-1
		}while(k==j);


		perm(j,k);
		mutated=true;
	}

	if(rnd->Rannyu()<prob_m){
		
		int j = 1 + int(rnd->Rannyu(0, N-2));		//intero tra 1 e N-2
		int k = 1 + int(rnd->Rannyu(0, N-2));		//intero tra 1 e N-2	
		
		int m = abs(k-j) + 1;	//numero di elementi da shiftare

		if(m>=N-1){		//m < N-1
			if(k>j)	k=N-2;
			else	j=N-2;
		}	

		int n = 1 + int(rnd->Rannyu(0, N-1-k));
		if(k<j)	n = 1 + int(rnd->Rannyu(0, N-1-j));


		shift(n,j,k); //k+n<N
		mutated=true;
	}

	if(rnd->Rannyu()<prob_m){

		int m_max=(N-1)/2;

		//int m = 1 + int(rnd->Rannyu(0, m_max));		//intero tra 1 e (N-1)/2
		
		int k = int(rnd->Rannyu(N-m_max, N));		//intero tra N-m_max e N-1 (nella seconda metà del vettore permutabile)
		int m = 1 + int(rnd->Rannyu(0, N-k));		//intero tra 1 e N-k (numero massimo di elementi del blocco)

		int j = 1 + int(rnd->Rannyu(0, k-m));		//intero tra 1 e k-m
		
		perm_m(m,j,k);
		mutated=true;
	}

	if(rnd->Rannyu()<prob_m){
		
		int j = 1 + int(rnd->Rannyu(0, N-1));	//intero tra 1 e N-1
		int k;

		do{ k = 1 + int(rnd->Rannyu(0, N-1));	//intero tra 1 e N-1
		}while(k==j);

		inversion(j,k);
		mutated=true;
	}


	if(mutated and !check())	reset(old);	//se check dà false aggiorno l'individuo con la copia non mutata

	return mutated;

}


double Individuo :: L1(vector<double> xc, vector<double> yc) const{

	if(xc.size()!=yc.size()){
		cerr<<"Coordinate città non valide (numero diverso di x e y)"<<endl;
		exit(-1);
	}

	int n_cities = xc.size();
	if(n_cities!=N){
		cerr<<"Numero di città incompatibile con la lunghezza del percorso"<<endl;
		exit(-2);
	}	
	
	double sum=0.;

	for(int i=0; i<N; i++){
		int pre, next;
		pre=path[i]-1;
		if(i==N-1)	next=path[0]-1;
		else		next=path[i+1]-1;

		sum+=sqrt( pow(xc[pre] - xc[next], 2) + pow(yc[pre] - yc[next], 2) );
	}

	return sum;
	

}


double Individuo :: L2(vector<double> xc, vector<double> yc) const{

	if(xc.size()!=yc.size()){
		cerr<<"Coordinate città non valide (numero diverso di x e y)"<<endl;
		exit(-1);
	}

	int n_cities = xc.size();
	if(n_cities!=N){
		cerr<<"Numero di città incompatibile con la lunghezza del percorso"<<endl;
		exit(-2);
	}

	double sum=0.;

	for(int i=0; i<N; i++){
		int pre, next;
		pre=path[i]-1;
		if(i==N-1)	next=path[0]-1;
		else		next=path[i+1]-1;

		sum+= pow(xc[pre] - xc[next], 2) + pow(yc[pre] - yc[next], 2);
	}

	return sum;
}


/********************************************** POPOLAZIONE ************************************************************/


Popolazione :: Popolazione(){

	InitRandom(&rnd);

	I=100;
	N=32;

	vector<Individuo> new_g(I);
	vector<Individuo> p(I);
	p[0]=Individuo(N);

	
	for(int i=1; i<I; i++){

		p[i]=Individuo(p[i-1]);		//copio l'individuo precedente

		int j = 1 + int(rnd.Rannyu(0, N-1));	//intero tra 1 e N-1
		int k = 1 + int(rnd.Rannyu(0, N-1));	//intero tra 1 e N-1

		p[i].perm(j, k);		//permuto due indici scelti a caso tra 1 e N-1

	}

	pop=p;
	new_gen=new_g;

	vector<double> x(N);
	vector<double> y(N);
	xc=x;
	yc=y;

	sorted=false;
	updated=false;


}




Popolazione :: Popolazione(int Is, int Ns){

	InitRandom(&rnd);

	I=Is;
	N=Ns;

	vector<Individuo> new_g;
	vector<Individuo> p(I);
	p[0]=Individuo(N);		//primo individuo [1, 2, ..., N]

	
	for(int i=1; i<I; i++){

		p[i]=Individuo(p[i-1]);		//copio l'individuo precedente

		int j = 1 + int(rnd.Rannyu(0, N-1));	//intero tra 1 e N-1
		int k = 1 + int(rnd.Rannyu(0, N-1));	//intero tra 1 e N-1

		p[i].perm(j, k);		//permuto due indici scelti a caso tra 1 e N-1

	}

	pop=p;
	new_gen=new_g;

	vector<double> x(N);
	vector<double> y(N);
	xc=x;
	yc=y;

	sorted=false;
	updated=false;

}


Popolazione :: Popolazione(Individuo ind){

	InitRandom(&rnd);

	I=100;
	N=ind.get_size();

	vector<Individuo> new_g(I);
	vector<Individuo> p(I);
	p[0]=ind;			//primo individuo passato

	
	for(int i=1; i<I; i++){

		p[i]=Individuo(p[i-1]);		//copio l'individuo precedente

		int j = 1 + int(rnd.Rannyu(0, N-1));	//intero tra 1 e N-1
		int k = 1 + int(rnd.Rannyu(0, N-1));	//intero tra 1 e N-1

		p[i].perm(j, k);		//permuto due indici scelti a caso tra 1 e N-1

	}

	pop=p;
	new_gen=new_g;

	vector<double> x(N);
	vector<double> y(N);
	xc=x;
	yc=y;

	sorted=false;
	updated=false;

}


Popolazione :: Popolazione(int Is, Individuo ind){

	InitRandom(&rnd);

	I=Is;
	N=ind.get_size();

	vector<Individuo> new_g(I);
	vector<Individuo> p(I);
	p[0]=ind;

	
	for(int i=1; i<I; i++){

		p[i]=Individuo(p[i-1]);		//copio l'individuo precedente

		int j = 1 + int(rnd.Rannyu(0, N-1));	//intero tra 1 e N-1
		int k = 1 + int(rnd.Rannyu(0, N-1));	//intero tra 1 e N-1

		p[i].perm(j, k);		//permuto due indici scelti a caso tra 1 e N-1

	}

	pop=p;
	new_gen=new_g;

	vector<double> x(N);
	vector<double> y(N);
	xc=x;
	yc=y;

	sorted=false;
	updated=false;


}


Popolazione :: Popolazione(const Popolazione &p_copy){


	I=p_copy.get_I();
	N=p_copy.get_ind(0).get_size();

	vector<Individuo> new_g(I);
	vector<Individuo> p(I);
	
	for(int i=0; i<I; i++)	p[i]=p_copy.get_ind(i);
	
	pop=p;
	new_gen=new_g;

	xc = p_copy.xc;
	yc = p_copy.yc;

	sorted=p_copy.sorted;
	updated=p_copy.updated;

	rnd = p_copy.rnd;

}


Popolazione :: ~Popolazione(){}


Popolazione& Popolazione::operator=(const Popolazione& p){

	I = p.I;
	N = p.N;
	
	if(!pop.empty()){
		pop.clear();
		pop.resize(0);
		pop.shrink_to_fit();
	}
	pop=p.pop;

	if(!new_gen.empty()){
		new_gen.clear();
		new_gen.resize(0);
		new_gen.shrink_to_fit();
	}
	new_gen=p.new_gen;

	if(!xc.empty()){
		xc.clear();
		xc.resize(0);
		xc.shrink_to_fit();
	}
	xc=p.xc;
	
	if(!yc.empty()){
		yc.clear();
		yc.resize(0);
		yc.shrink_to_fit();
	}
	yc=p.yc;

	sorted = p.sorted;
	updated = p.updated;
	rnd = p.rnd;

	return *this;

}




/****************************************************** METODI **********************************************************/


void Popolazione :: set_cities(vector<double> x, vector<double> y){


	if(x.size()!=y.size()){
		cerr<<"Coordinate città non valide (numero diverso di x e y)"<<endl;
		exit(-1);
	}

	int n_cities = x.size();
	if(n_cities!=N){
		cerr<<"Numero di città incompatibile con la lunghezza del percorso"<<endl;
		exit(-2);
	}

	xc=x;
	yc=y;


}


Individuo Popolazione :: get_ind(int i) const{

	if(i<0){

		cerr<<"L'indice richiesto è negativo. Returning i ="<<0<<endl;
		return pop[0];

	}	

	if(i>=I){

		cerr<<"L'indice richiesto non è nella popolazione. Returning i ="<<I-1<<endl;
		return pop[I-1];

	}

	else return pop[i];

}



bool Popolazione :: check_ind(int i){

	Individuo ind;
	ind=pop[i];

	bool ok=ind.check();	//la copia di pop[i] nel caso viene modificata


	if(!ok){		//se l'individuo non va bene lo cambio con

		sorted=false;
		updated=false;

	}

	pop[i].reset(ind);	//resetto l'individuo in ogni caso, perché anche se ok=true l'individuo potrebbe essere stato ruotato, ma la fitness è uguale
	

	return ok;

}

void Popolazione :: printPop() const{

	for(int j=0; j<I; j++)		pop[j].printInd();


}

void Popolazione :: printNG() const{

	int size=new_gen.size();
	
	for(int j=0; j<size; j++){
		new_gen[j].printInd();
	}


}


void Popolazione :: update_fitness(){

	if(updated)	return;

	for(int i=0; i<I; i++)	pop[i].set_fitness(xc, yc);

	updated=true;
	sorted=false;


}

void Popolazione :: update_fitness(int i){


	pop[i].set_fitness(xc, yc);

	sorted=false;

}

void Popolazione :: update_fitness_ng(){

	for(int i=0; i<I; i++)	new_gen[i].set_fitness(xc, yc);


}

void Popolazione :: update_fitness_ng(int i){

	new_gen[i].set_fitness(xc, yc);

}



void Popolazione :: sort_pop(){

	update_fitness();
	
	//sort(pop.begin(), pop.end(), [](const Individuo & l, const Individuo & r){ return l.get_fitness() < r.get_fitness(); } );
	sort(pop.begin(), pop.end());

	sorted=true;

}

void Popolazione::quicksort(int primo, int ultimo){


	if(primo>ultimo or ultimo>=I) return;

	update_fitness();
	
	if(ultimo-primo<=1){
	
		if(get_fitness(ultimo)<get_fitness(primo))
			scambia(primo, ultimo);

		return;
		

	}
	
	double pivot=get_fitness((ultimo+primo)/2);
	int basso=primo;
	int alto=ultimo;
	
	while(basso<alto){
	
		while(get_fitness(basso)<pivot) basso++;
		while(get_fitness(alto)>pivot) alto--;
		
		if(basso<alto){
		scambia(basso, alto);
		basso++;
		}
	
	}


	quicksort(primo, basso-1);
	quicksort(basso, ultimo);
	
	sorted=true;

}


void Popolazione :: scambia(int i, int j){

	if(i<0 or i>=I or j<0 or j>=I){
		cout<<"Gli indici "<<i<<", "<<j<<" sono fuori dal range consentito, permutazione non eseguita"<<endl;
		return;
	}

	Individuo a(get_ind(i));
	pop[i]=pop[j];
	pop[j]=a;

	sorted=false;


}

int Popolazione :: selection_unif(){

	if(!sorted)	quicksort(0, I-1);


	double r=rnd.Rannyu();
	int p=5;

	return int(I* pow(r, p));		

}


int Popolazione :: selection(){

	if(!sorted)	quicksort(0, I-1);

	vector<double> prob(I);

	double prob_tot=0;

	for(int i=0; i<I; i++){
		
		if(i==0)	prob[i]=1;
		else{
			double delta = get_fitness(i) - get_fitness(0);
			prob[i]=pow(delta+1, -2);		//più la fitness si discosta dal primo, più bassa è la probabilità, delta +1 per evitare valori >1
			
		}

		prob_tot+=prob[i];

	}


	double r=rnd.Rannyu(0, prob_tot);

	double prob_part=0;		//probabilità parziale = somma delle probabilità testate fino a quel punto, dà il limite inferiore di ricerca
	int index=0;

	for(int i=0; i<I; i++){

		
		if(prob_part<r and r<prob_part + prob[i]){
			index=i;
			break;
		}

		prob_part+=prob[i];

	}

	return index;		

}



void Popolazione :: mutate_ind(int i, double prob_m){

	Individuo ind;
	ind=pop[i];

	bool mutated=ind.mutate(prob_m, &rnd);

	if(mutated){	//se mutated dà false lascio l'individuo vecchio (la funzione mutate controlla anche la validità dell'individuo mutato)
		pop[i].reset(ind);	//aggiorno l'individuo con la copia mutata
		update_fitness(i);	//aggiorno la fitness dell'individuo mutato
		sorted=false;
	}


}


void Popolazione :: mutate_all(double prob_m){

	for(int i=0; i<I; i++)	mutate_ind(i, prob_m);
	
	updated=true;

}


void Popolazione :: crossover (int i, int j, double prob_c){

	bool crossed=false;

	vector<int> madre(get_ind(i).get_vec());
	vector<int> padre(get_ind(j).get_vec());

	if(rnd.Rannyu()<prob_c){

		int n = 1 + rnd.Rannyu(0, N-1);	//indice tra 1 e N-1
		

		vector<int> cross_m(N-n);
		vector<int> cross_p(N-n);

		for(int k=0; k<N-n; k++){
			cross_m[k]=madre[n+k];
			cross_p[k]=padre[n+k];
		}

		for(int k=1; k<N; k++){		//scorro i due genitori per trovare gli elementi che sono nella parte da incrociare

			for(int l=0; l<N-n; l++){

				if(padre[k]==cross_m[l]){	//ordino la parte da incorciare della madre in modo che i suoi elementi compaiano come nel padre
					
					cross_m.erase(cross_m.begin() + l);
					cross_m.insert(cross_m.end(), padre[k]);
				}

				if(madre[k]==cross_p[l]){	//ordino la parte da incorciare del padre in modo che i suoi elementi compaiano come nella madre
					
					cross_p.erase(cross_p.begin() + l);
					cross_p.insert(cross_p.end(), madre[k]);
				}

			}

		}

	
		for(int k=0; k<N-n; k++){
			madre[n+k]=cross_m[k];
			padre[n+k]=cross_p[k];
		}

		crossed=true;

	}

	Individuo figlio1(madre);
	Individuo figlio2(padre);

	if(crossed){

		if(!figlio1.check())	figlio1=get_ind(i);
		if(!figlio2.check())	figlio2=get_ind(j);


	}

	figlio1.set_fitness(xc, yc);
	figlio2.set_fitness(xc, yc);

	new_gen.push_back(figlio1);
	new_gen.push_back(figlio2);



}


void Popolazione :: update_gen(){

	int size_ng=new_gen.size();


	if(size_ng>I){

		cerr<<"Errore: nuova generazione più numerosa della precedente"<<endl;
		exit(-3);
	
	}

	
	if(size_ng<I){

		update_fitness();
		if(!sorted)	quicksort(0, I-1);
		for(int i=0; i<I-size_ng; i++)	new_gen.push_back(pop[i]);	//per gli elementi mancanti prendo i più fit dalla vecchia generazione

	}

	pop=new_gen;
	update_fitness();

	new_gen.clear();
	new_gen.resize(0);
	new_gen.shrink_to_fit();

}


void Popolazione :: RS(double prob_m){


	update_fitness();
	mutate_all(prob_m);



}


void Popolazione :: GA(double prob_m, double prob_c){


	update_fitness();

	for(int i=0; i<I/2; i++){

		int j=selection();
		int k;

		do{ k=selection(); }while(k==j);

		crossover(j, k, prob_c);
		//cout<<i<<" "<<new_gen.size()<<endl;			
		
	}

	update_gen();

	mutate_all(prob_m);



}


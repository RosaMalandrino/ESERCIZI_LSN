#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "random.h"
#include "genetic.h"

using namespace std;


double error(double av[], double av2[], int n);


int main (int argc, char *argv[]){



/********************************INIZIALIZZAZIONE GENERATORE DI NUMERI CASUALI**********************************/
	
	Random rnd;
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
            			rnd.SetRandom(seed,p1,p2);
         		}
      		}
      	input.close();
   	} else cerr << "PROBLEM: Unable to open seed.in" << endl;



/****************************************ESERCIZIO 09.1**************************************************************/


	int N=32;	//numero di città

	vector<double> xs(N);
	vector<double> ys(N);
	

	bool circle=false;	//true per il cerchio, false per il quadrato

	string suffix;
	if(circle)	suffix = "_circle.dat";
	else		suffix = "_square.dat";


	ifstream cities("cities" + suffix);

	for(int i=0; i<N; i++)	cities>>xs[i]>>ys[i];

	cities.close();


	Individuo path(N);		//parto dal percorso [1, ..., N]
	path.set_fitness(xs, ys);
		
	ofstream fitness, LofT, path_opt;

	fitness.open("Ls" + suffix);
	LofT.open("LofT" + suffix);
	path_opt.open("optimal_path" + suffix);


	double p_m=0.7;

	int n=10000;		//suddivisioni della temperatura
	int n_min=100;		//numero minimo  di iterazioni per temperatura
	int n_max=1000;		//numero massimo di iterazioni per temperatura

	double Tmax=2.;
	double dT=Tmax/n;
	double T=Tmax;
	
	


	for(int j=0; j<n; j++){

		T=Tmax-j*dT;
		double beta=1/T;
		
		int n_acc=0;
		int n_tot=0;
		double L_old=0., L_new=0.;
		double L=path.get_fitness();
		double L_pre=L;
		int n_same=0;


		for(int i=0; i<n_max; i++){

			if(i%n_min==0)	n_same=0;		//azzero il conteggio delle ripetizioni ogni n_min volte
		
			L_old = path.get_fitness();
			Individuo path_new(path);		//copio il vecchio percorso
			path_new.mutate(p_m, &rnd);		//provo a mutarlo
			path_new.set_fitness(xs, ys);		//calcolo la fitness del nuovo percorso
			L_new = path_new.get_fitness();

			double p=exp(-beta*(L_new - L_old));

			if(rnd.Rannyu()<=p){

				path=path_new;
				L=path.get_fitness(); 
				n_acc++;

			}

			n_tot++;


			if(L==L_pre)	n_same++;	//conto quante volte il nuovo path è uguale a quello dell'iterazione precedente
			L_pre=L;			//aggiorno il valore della fitness precedente per la prossima iterazione

		
			fitness<<L<<endl;

			lasti=i;
			if(i%n_min==n_min-1 and n_same<n_min)	break;		//se in n_min iterazioni ci sono stati cambiamenti esco e abbasso la temperatura
			

		}

		LofT<<T<<" "<<L<<endl;		//salvo il valore della fitness in funzione della temperatura
		cout<<"Temperatura: "<<T<<endl;
		cout<<"Rate di accettazione: "<<((double) n_acc)/n_tot<<endl<<endl;
		
	

	
	}

	path.printInd();

	for(int i=0; i<path.get_size(); i++)	path_opt<<path.get_comp(i)<<endl;
		

	fitness.close();
	LofT.close();
	path_opt.close();


	return 0;

}



double error(double av[], double av2[], int n){


	if (n == 0){
		return 0;
	}
	else{
		 return sqrt((av2[n]- (av[n])*(av[n]))/n);
	}		


}

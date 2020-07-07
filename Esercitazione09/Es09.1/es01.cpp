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
	int I=100;	//numero di individui

	vector<double> xs(N);
	vector<double> ys(N);
	

	bool circle=false;	//true per il cerchio, false per il quadrato

	string suffix;
	if(circle)	suffix = "_circle.dat";
	else		suffix = "_square.dat";


	ofstream cities, best, ave, path;

	cities.open("cities" + suffix);

	if(circle){

		for(int i=0; i<N; i++){

			double x=0., y=0., r=0.;

			do{					//generazione di punti casuali all'interno della circonferenza di raggio 1 centrato in (0, 0)
				x=rnd.Rannyu(-1, 1);
				y=rnd.Rannyu(-1, 1);
				r=sqrt(pow(x,2)+pow(y,2));

			}while(r>1);

			xs[i]=x/r;
			ys[i]=y/r;

			cities<<xs[i]<<" "<<ys[i]<<endl;	
		}

	}		



	if(!circle){

		for(int i=0; i<N; i++){
								//generazione di punti casuali all'interno della quadrato di lato 2 centrato in (0, 0)
			xs[i]=rnd.Rannyu(-1, 1);
			ys[i]=rnd.Rannyu(-1, 1);

			cities<<xs[i]<<" "<<ys[i]<<endl;
		}

	}

	cities.close();


	Popolazione p(I, N);
	p.set_cities(xs, ys);

		

	best.open("best" + suffix);
	ave.open("ave" + suffix);
	path.open("optimal_path" + suffix);


	//p.update_fitness();
		

	if(!p.IsSorted())	p.quicksort(0, p.get_I()-1);

	double lowest=p.get_fitness(0);
	
	Individuo ind_opt(N);


	int n_iter=500;

	if(circle)	n_iter=1000;

	for(int j=0; j<n_iter; j++){

		p.GA(0.05, 0.75);	//step dell'algoritmo genetico, con probabilità di mutazione del 5% e di crossover del 75%
		//p.RS(0.7);


		if(!p.IsSorted())	p.quicksort(0, p.get_I()-1);

		double L_b=p.get_fitness(0);	//fitness del miglior individuo

		best<<L_b<<endl;

			

		double sum=0.;

		for(int i=0; i<p.get_I()/2; i++)	sum+=p.get_fitness(i);
	
		ave<<2*sum/p.get_I()<<endl;	//fitness mediata sulla metà migliore della popolazione
			
		if(L_b<lowest){
			lowest=L_b;
			ind_opt=p.get_ind(0);
		}

	}


	for(int i=0; i<ind_opt.get_size(); i++)	path<<ind_opt.get_comp(i)<<endl;


	cout<<"Best path in final population"<<endl;
	p.get_ind(0).printInd();

	cout<<"Optimal path"<<endl;
	ind_opt.printInd();

	best.close();
	ave.close();
	path.close();


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

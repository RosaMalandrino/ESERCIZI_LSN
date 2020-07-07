#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "random.h"
#include "genetic.h"

#include "mpi.h"

using namespace std;


double error(double av[], double av2[], int n);


int main (int argc, char *argv[]){

	
	int size, rank;
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);


/********************************INIZIALIZZAZIONE GENERATORE DI NUMERI CASUALI**********************************/
	
	Random rnd;
	int seed[4];
	int p1, p2;
	ifstream Primes("Primes");
	if (Primes.is_open()){
		for(int i=0; i<rank+1; i++)	Primes >> p1 >> p2 ;
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



/****************************************ESERCIZIO 10.2**************************************************************/


	int N=32;	//numero di città
	int I=100;	//numero di individui

	vector<double> xs(N);
	vector<double> ys(N);
	

	ifstream cities("cities_square.dat");
	
	for(int i=0; i<N; i++)		cities>>xs[i]>>ys[i];


	cities.close();

	

	Popolazione p(I, N);
	p.set_cities(xs, ys);
	p.set_random(rnd);	
		
	
	if(!p.IsSorted())	p.quicksort(0, p.get_I()-1);

	Individuo best_ind=p.get_ind(0);

	double lowest=p.get_fitness(0);

	Individuo ind_opt(N);


	int n_iter=500;
	int N_gen=50;


	ofstream best, ave;

	best.open("best_square_cont" + to_string(rank+1) + "_"+ to_string(N_gen) + ".dat");
	ave.open("ave_square_cont" + to_string(rank+1) + "_"+ to_string(N_gen) + ".dat");
	

	for(int j=0; j<n_iter; j++){

		p.GA(0.05, 0.75);	//step dell'algoritmo genetico, con probabilità di mutazione del 5% e di crossover del 75%

		if(!p.IsSorted())	p.quicksort(0, p.get_I()-1);

		best_ind=p.get_ind(0);

		double L_b=p.get_fitness(0);	//fitness del miglior individuo
		best<<L_b<<endl;

			

		double sum=0.;

		for(int i=0; i<p.get_I()/2; i++)	sum+=p.get_fitness(i);		

		ave<<2*sum/p.get_I()<<endl;					//fitness mediata sulla metà migliore della popolazione
			
		if(L_b<lowest){
			lowest=L_b;
			ind_opt=p.get_ind(0);
		}

	//MIGRAZIONE

		if(j==n_iter-1)	break;
	
		if(j%N_gen==N_gen-1){


			int *oldpath=new int[N];
			int *newpath=new int[N];

			for(int k=0; k<N; k++)	oldpath[k]=best_ind.get_comp(k);

			int cont1=-1, cont2=-1, cont=-1;		//in questo modo nessun rank può essere uguale agli indici di continenti

			
			if(rank==0){				//estraggo casualmente il rank da accoppiare a 0 (cont1=0, cont2=cont)

				cont1=0;

				do{
					cont=int(size*rnd.Rannyu());
				}while(cont==cont1);

				cont2=cont;

			}

			MPI_Bcast(&cont, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);	//trasmetto agli altri nodi il rank estratto


			if(rank==cont){		//dico al continente associato a 0 chi è cont1 e cont2

				cont1=0;
				cont2=cont;

			}

			if(rank!=cont and rank!=0){	//dico agli altri due continenti chi è cont1 e chi cont2

				switch(cont){

					case(1):
						cont1=2;
						cont2=3;
						break;

					case(2):
						cont1=1;
						cont2=3;
						break;

					case(3):
						cont1=1;
						cont2=2;
						break;
				}

			}


			MPI_Status stat1, stat2;
			MPI_Request req;

			int itag=1;
			int itag2=2;

			if(rank==cont1){

				MPI_Isend(&oldpath[0], N, MPI_INTEGER, cont2, itag, MPI_COMM_WORLD, &req);
				MPI_Recv(&newpath[0], N, MPI_INTEGER, cont2, itag2, MPI_COMM_WORLD, &stat1);

			}

			else if(rank==cont2){
				
				MPI_Send(&oldpath[0], N, MPI_INTEGER, cont1, itag2, MPI_COMM_WORLD);
				MPI_Recv(&newpath[0], N, MPI_INTEGER, cont1, itag, MPI_COMM_WORLD, &stat2);
				
			}
			
			
			Individuo ind(N);

			for(int k=0; k<N; k++)	ind.set_comp(k, newpath[k]);
			ind.set_fitness(xs, ys);

			p.set_ind(0, ind);

			delete [] oldpath;
			delete [] newpath;

		}


		
	}





/*	
	cout<<"Best path in final population"<<endl;
	p.get_ind(0).printInd();

	cout<<"Optimal path"<<endl;
	ind_opt.printInd();
*/

	cout<<"Continente "<<rank<<":"<<endl;
	ind_opt.printInd();
	cout<<endl;

	double *fitness=new double[size];

	for(int i=0; i<size; i++){

		if(i==rank)	fitness[i]=ind_opt.get_fitness();
	
		MPI_Bcast(&fitness[i], 1, MPI_REAL8, i, MPI_COMM_WORLD);	
	}

	double best_fit=fitness[0];
	int best_cont=0;

	for(int i=0; i<size; i++){	
		if(fitness[i]<best_fit){
			best_fit=fitness[i];
			best_cont=i;
		}

	}

	//if(rank==0)	cout<<best_cont<<endl;


	ofstream path, best_path;

	best_path.open("optimal_path_square_" + to_string(N_gen) + ".dat");
	path.open("optimal_path_square_cont"+ to_string(rank+1) +"_" + to_string(N_gen) + ".dat");

	for(int i=0; i<ind_opt.get_size(); i++){
		path<<ind_opt.get_comp(i)<<endl;
		if(rank==best_cont)	best_path<<ind_opt.get_comp(i)<<endl;
	}

	best.close();
	ave.close();
	path.close();
	best_path.close();



	MPI_Finalize();

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

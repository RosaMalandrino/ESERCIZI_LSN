#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

using namespace std;


double error(double av[], double av2[], int n);
double gauss(double x, double mu, double sigma){ return exp(-0.5*pow((x-mu)/sigma,2)); }
double psi(double x, double mu, double sigma){ return gauss(x, mu, sigma) + gauss(x, -mu, sigma); }
double V(double x){ return pow(x,4) -2.5*pow(x,2);}
double H(double x, double mu, double sigma);


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



/****************************************VARIATIONAL MONTE CARLO**************************************************************/


	//NOTA: in fondo al codice vengono salvati i valori finale sul file ottimizzazione.dat. Commentare quando non si sta ottimizzando!


	int M=100000;			//numero di step (tenere <=100000 per ottimizzazione con simulated annealing)
	int N=100;			//blocchi
	
	

	//Valori di mu e sigma da linea di comando


	if(argc!=3){
		cout << "Input non valido." << endl;
		cout << "Comando corretto: " << argv[0] << " [mu] [sigma] " << endl;
		return -1;
	}

	double mu=atof(argv[1]);
	double sigma=atof(argv[2]);


	

	int n_eq=100;	//numero di step per equilibrare, è sufficiente anche per i punti di partenza molto lontani
	int dim=n_eq+M;

	double *xs=new double[2*dim];
	double *Hs=new double[M];

	//posizione iniziale

	double x0=mu+sigma;
	xs[0]=x0;
	xs[1]=-x0;

	int n_acc=0;		//passi accettati
	int n_tot=0;		//passi totali	

	double step_test=2*sigma;
	

	//if(sigma/mu <0.5)	step_test=2*mu;		//i picchi sono troppo stretti e distanti, serve un passo molto largo per campionarli entrambi
	

	//assegnazione del step per tenere l'accettazione tra il 40 e il 60 %

	if(sigma/mu > 0.55 and sigma/mu < 1.0)	step_test=4*sigma;
	if(sigma/mu > 1.0 and sigma/mu < 2.5)	step_test=3*sigma;		
 

	ofstream camp;
	camp.open("campionamento.dat");

	

	//sfrutto la simmetria della funzione, per ogni punto x campionato aggiungo anche -x


	for(int i=1; i<dim; i++){

			double xold=xs[2*i-2];
			
			double xtest= xold + step_test*rnd.Rannyu(-1, 1);		//step uniforme
				
			double psi_test, psi_old, ptest, pold, A;

			psi_test = psi(xtest, mu, sigma);		//funzione psi_T valutata nel punto di prova
			psi_old = psi(xold, mu, sigma);			//funzione psi_T valutata nel punto vecchio

			ptest=pow(psi_test, 2);				//probabilità del punto di prova
			pold =pow(psi_old , 2);				//probabilità del punto vecchio

			A = ptest/pold;
			
			double rand=rnd.Rannyu();

			if(rand<=A){		//accetto
				xs[2*i]=xtest;
				n_acc++;
				
			}
			else{			//rigetto
				xs[2*i]=xold;
			}

			n_tot++;

			xs[2*i+1]=-xs[2*i];


			if(i>=n_eq){
				Hs[i-n_eq] = H(xs[2*i], mu, sigma);		//valuto H solo dopo l'equilibrazione, lo faccio solo in xs, perché in -xs è uguale
				camp<<xs[2*i]<<endl<<xs[2*i+1]<<endl;
			}

		}

		cout<<"Rate di accettazione: "<<((double) n_acc)/n_tot<<endl;


	


	camp.close();


	ofstream fileout, opt;
	fileout.open("risultati.dat");
	opt.open("ottimizzazione.dat", ios::app);

	double H_ave[N];		//stima del valor medio <r>
	double H2_ave[N];		//stima del valor medio <r>^2

	double sum_prog[N]={};	
	double sum2_prog[N]={};	
	double err_prog[N];

	
	int L=int(M/N);

	for(int i=0; i<N; i++){
		
		double sum=0;
	
		for(int j=0; j<L; j++){

			int k = j + i*L;

			sum+=Hs[k];
				
		}

			
		H_ave[i]=sum/L;			//media dei valori di H psi/psi sui punti campionati secondo psi^2
		H2_ave[i]=pow(H_ave[i], 2);

		
		for(int j=0; j<i+1; j++){

			sum_prog[i]+=H_ave[j];	
			sum2_prog[i]+=H2_ave[j];
					
		}

		sum_prog[i]/=(i+1);
		sum2_prog[i]/=(i+1);

			
		err_prog[i] = error(sum_prog,sum2_prog,i);

		fileout<<(i+1)*int(M/N)<<" "<<sum_prog[i]<<" "<<err_prog[i]<<endl;	//step Monte Carlo
		//fileout<<(i+1)*L<<" "<<sum_prog[i]<<" "<<err_prog[i]<<endl;		//

	}

	fileout.close();

	//opt<<sum_prog[N-1]<<" "<<mu<<" "<<sigma<<endl;	//commentare quando non si sta ottimizzando

	opt.close();


	delete []xs;
	delete []Hs;



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




double H(double x, double mu, double sigma){


	double g_meno = gauss(x, mu, sigma);		//gaussiana con (x-mu)
	double g_piu = gauss(x, -mu, sigma);		//gaussiana con (x+mu)


	double H_psi = 0.5*pow(sigma, -2) * ( (1 - pow((x-mu)/sigma , 2))*g_meno + (1 - pow((x+mu)/sigma , 2))*g_piu ); // + V(x)*psi(x, mu, sigma)

	return H_psi/psi(x, mu, sigma) + V(x);
}

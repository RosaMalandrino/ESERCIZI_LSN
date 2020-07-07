#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

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



/****************************************ESERCIZIO 01.3**************************************************************/

	
	int M=1000000;		//numero di lanci
	int N=100;		// Numero di blocchi (esperimenti)
	int K=int(M/N);		// Numero di lanci per blocco (M deve essere multiplo di N)

	double L=0.3;		//lunghezza ago
	double d=1.;		//distanza righe

	double pi[N];		// Vettore delle stime di pi in ciascun blocco
	double pi2[N];		// Vettore delle stime di pi^2 in ciascun blocco al quadrato
	double sum_prog[N];	// Vettore delle somme delle stime di pi al progredire degli esperimenti
	double sum2_prog[N];	// Vettore delle somme delle stime di pi^2 al progredire degli esperimenti
	double err_prog[N];	// Vettore degli errori al progredire degli esperimenti



	double *x1=new double[M];	//vettore delle coordinate x del primo estremo dell'ago
	double *x2=new double[M];	//vettore delle coordinate x del secondo estremo dell'ago



	for(int i=0; i<M; i++){

		double xc=rnd.Rannyu(-d/2, d/2);	//coordinata x del centro dell'ago nell'intervallo [-d/2, d/2]

		double x=0;
		double y=0;
		double r=0;
		
		do{					//generazione di punti casuali all'interno della circonferenza di raggio L/2
			x=rnd.Rannyu(-L/2, L/2);
			y=rnd.Rannyu(-L/2, L/2);
			r=sqrt(pow(x,2)+pow(y,2));

		}while(r>0.5*L);
		
		double dx=0.5*L*x/r;				//componente x della proiezione di questi punti sulla circonferenza (distribuiti uniformemente)


		//double theta=rnd.Rannyu(0., M_PI);		//calcolo alternativo generando l'angolo teta tra l'ago e l'asse x nell'intervallo [0, pi greco] 
		//dx=0.5*L*cos(theta);			


		x1[i]=xc + dx;
		x2[i]=xc - dx;



	}

	

	


	for(int i=0; i<N; i++){
	
		int Nhit=0;
		int Ntot=0;

		for(int j=0; j<K; j++){

			int k = j + i*K;
				
 			if(x1[k]*x2[k]<0.)	Nhit++;		//se x1 e x2 hanno segno opposto l'ago attraversa la linea in x=0

			Ntot++;
		}

		pi[i] = 2*L*Ntot/(Nhit*d);
		pi2[i] = pow(pi[i], 2);
		
	}



	ofstream fileout("risultati_01.3.dat");

	for(int i=0; i<N; i++){
	
		for(int j=0; j<i+1; j++){

			sum_prog[i]+=pi[j];	//somma progressiva delle stime di pi fino all'i-esimo esperimento
			sum2_prog[i]+=pi2[j];	//somma progressiva delle stime di pi^2 fino all'i-esimo esperimento
		}

		sum_prog[i]/=(i+1);		//media delle stime di pi fino all'i-esimo esperimento
		sum2_prog[i]/=(i+1);		//media delle stime di pi^2 fino all'i-esimo esperimento

		err_prog[i] = error(sum_prog,sum2_prog,i);	//errore nella stima della media all'i-esimo esperimento

	
		fileout<<(i+1)*K<<" "<<sum_prog[i]<<" "<<err_prog[i]<<endl;

	}

	fileout.close();

	delete [] x1;
	delete [] x2;


	
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

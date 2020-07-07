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
//p1=279;
//p2=39534983;
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



/****************************************ESERCIZIO 03.1**************************************************************/


	int M=1000000;		//numero di simulazioni
	int N=100;		//numero di blocchi
	int L=int(M/N);

	int n=100;		//suddivisioni temporali


	double S0=100.;		//prezzo iniziale
	double T=1.;		//delivery time
	double K=100.;		//strike price
	double r=0.1;		//tasso d'interesse
	double sigma=0.25;	//volatilità

	double dt=T/n;


	double C[N];		//stime del costo della call per ogni blocco   (metodo diretto)
	double C2[N];		//stime^2 del costo della call per ogni blocco (metodo diretto)

	double P[N];		//stime del costo della put per ogni blocco    (metodo diretto)
	double P2[N];		//stime^2 del costo della put per ogni blocco  (metodo diretto)

	double sum_prog_C[N];
	double sum2_prog_C[N];
	double err_prog_C[N];

	double sum_prog_P[N];
	double sum2_prog_P[N];
	double err_prog_P[N];



	double C_dis[N];	//stime del costo della call per ogni blocco   (con discretizzazione)
	double C2_dis[N];	//stime^2 del costo della call per ogni blocco (con discretizzazione)

	double P_dis[N];	//stime del costo della put per ogni blocco    (con discretizzazione)
	double P2_dis[N];	//stime^2 del costo della put per ogni blocco  (con discretizzazione)

	double sum_prog_C_dis[N];
	double sum2_prog_C_dis[N];
	double err_prog_C_dis[N];

	double sum_prog_P_dis[N];
	double sum2_prog_P_dis[N];
	double err_prog_P_dis[N];

	

	
	ofstream fileout("risultati_03.1.dat");

	for(int i=0; i<N; i++){
	
		double sum_C=0;		//accumulazione del costo della call
		double sum_P=0;		//accumulazione del costo della put

		double sum_C_dis=0;	//accumulazione del costo della call con discretizzazione
		double sum_P_dis=0;	//accumulazione del costo della put con discretizzazione


		for(int j=0; j<L; j++){


			double S=S0*exp( (r-0.5*pow(sigma, 2))*T + sigma*rnd.Gauss(0,sqrt(T)) );	//calcolo del prezzo finale S(T) con il metodo diretto


			double g_C=0;		//guadagno nel caso della call
			double g_P=0;		//guadagno nel caso della put

			if((S-K)>0)	g_C=S-K;
			else		g_P = K-S;

			sum_C+=exp(-r*T)*g_C;
			sum_P+=exp(-r*T)*g_P;


			// Con discretizzazione 

			double S_pre=S0;	//variabile per salvare il valore di S al tempo precedente, inizialmente =S0

			for(int k=0; k<n; k++){

				S_pre*=	exp( (r-0.5*pow(sigma, 2))*dt + sigma*sqrt(dt)*rnd.Gauss(0,1) );	//evoluzione di S dal punto precedente

			}

			S=S_pre;		//S al tempo finale T


			if((S-K)>0){
				g_C=S-K;
				g_P=0;
			}

			else{
				g_P = K-S;
				g_C=0;
			}

			sum_C_dis+=exp(-r*T)*g_C;
			sum_P_dis+=exp(-r*T)*g_P;

			
		}

		C[i]=sum_C/L;			
		C2[i]=pow(C[i], 2);

		P[i]=sum_P/L;			
		P2[i]=pow(P[i], 2);


		C_dis[i]=sum_C_dis/L;			
		C2_dis[i]=pow(C_dis[i], 2);

		P_dis[i]=sum_P_dis/L;			
		P2_dis[i]=pow(P_dis[i], 2);


		//calcolo di medie ed errori progressivi

		for(int j=0; j<i+1; j++){

			sum_prog_C[i]+=C[j];	
			sum2_prog_C[i]+=C2[j];

			sum_prog_P[i]+=P[j];	
			sum2_prog_P[i]+=P2[j];


			sum_prog_C_dis[i]+=C_dis[j];	
			sum2_prog_C_dis[i]+=C2_dis[j];

			sum_prog_P_dis[i]+=P_dis[j];	
			sum2_prog_P_dis[i]+=P2_dis[j];
	
		}

		sum_prog_C[i]/=(i+1);
		sum2_prog_C[i]/=(i+1);
		err_prog_C[i] = error(sum_prog_C,sum2_prog_C,i);

		sum_prog_P[i]/=(i+1);
		sum2_prog_P[i]/=(i+1);
		err_prog_P[i] = error(sum_prog_P,sum2_prog_P,i);


		sum_prog_C_dis[i]/=(i+1);
		sum2_prog_C_dis[i]/=(i+1);
		err_prog_C_dis[i] = error(sum_prog_C_dis,sum2_prog_C_dis,i);

		sum_prog_P_dis[i]/=(i+1);
		sum2_prog_P_dis[i]/=(i+1);
		err_prog_P_dis[i] = error(sum_prog_P_dis,sum2_prog_P_dis,i);


		fileout<<(i+1)*L<<" "<<sum_prog_C[i]<<" "<<err_prog_C[i]<<" "<<sum_prog_P[i]<<" "<<err_prog_P[i]<<" "
			 <<sum_prog_C_dis[i]<<" "<<err_prog_C_dis[i]<<" "<<sum_prog_P_dis[i]<<" "<<err_prog_P_dis[i]<<endl;

	
	}

	fileout.close();



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

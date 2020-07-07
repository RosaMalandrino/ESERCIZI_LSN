#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

using namespace std;


double error(double [], double [], int);


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



/****************************************ESERCIZIO 01.1.1**************************************************************/


	int M=1000000;		// Numero di lanci totali
	int N=100;		// Numero di blocchi (esperimenti)
	int L=int(M/N);		// Numero di lanci per blocco (M deve essere multiplo di N)
	
	double *rand=new double[M];	// Vettore di numeri casuali distribuiti uniformente tra 0 e 1
	double ave[N];			// Vettore delle medie dei lanci in ciascun blocco
	double ave2[N];			// Vettore delle medie dei lanci in ciascun blocco al quadrato
	double sum_prog[N];		// Vettore delle somme delle medie al progredire degli esperimenti
	double sum2_prog[N];		// Vettore delle somme delle medie al quadrato al progredire degli esperimenti
	double err_prog[N];		// Vettore degli errori al progredire degli esperimenti

	

	for(int i=0; i<M; i++){

		rand[i] = rnd.Rannyu();
		
	}

	rnd.SaveSeed();




	for(int i=0; i<N; i++){
	
		double sum=0;

		for(int j=0; j<L; j++){

			int k = j + i*L;
			sum += rand[k];
		}

		ave[i] = sum/L;
		ave2[i] = ave[i]*ave[i];
		
	}



	ofstream fileout1("risultati_01.1.1.dat");


	for(int i=0; i<N; i++){
	
		for(int j=0; j<i+1; j++){

			sum_prog[i]+=ave[j];	//somma progressiva delle medie fino all'i-esimo esperimento
			sum2_prog[i]+=ave2[j];	//somma progressiva delle medie^2 fino all'i-esimo esperimento
		}

		sum_prog[i]/=(i+1);		//media delle medie fino all'i-esimo esperimento
		sum2_prog[i]/=(i+1);		//media delle medie^2 fino all'i-esimo esperimento

		err_prog[i] = error(sum_prog,sum2_prog,i);	//errore nella stima della media all'i-esimo esperimento

	
		fileout1<<(i+1)*L<<" "<<sum_prog[i]<<" "<<err_prog[i]<<endl;

	}

	fileout1.close();



/****************************************ESERCIZIO 01.1.2**************************************************************/




	for(int i=0; i<N; i++){
	
		double sum=0;

		for(int j=0; j<L; j++){

			int k = j + i*L;
			sum += pow((rand[k]-0.5), 2);	//scarti quadratici dalla media
		}

		ave[i] = sum/L;			//sovrascrivo il vettore usato precedentemente per calcolare le stime della varianza
		ave2[i] = ave[i]*ave[i];	// 	""	""		""	""		""	""	 	""   al quadrato
		
	}



	ofstream fileout2("risultati_01.1.2.dat");	

	for(int i=0; i<N; i++){

		sum_prog[i]=0;		//azzero il vettore usato precedentemente
		sum2_prog[i]=0;		// 	""	""	""	""
	
		for(int j=0; j<i+1; j++){

			sum_prog[i]+=ave[j];	//somma progressiva delle varianze fino all'i-esimo esperimento
			sum2_prog[i]+=ave2[j];	//somma progressiva delle varianze^2 fino all'i-esimo esperimento
		}

		sum_prog[i]/=(i+1);		//media delle varianze fino all'i-esimo esperimento
		sum2_prog[i]/=(i+1);		//media delle varianze^2 fino all'i-esimo esperimento

		err_prog[i] = error(sum_prog,sum2_prog,i);	//errore nella stima della varianza all'i-esimo esperimento

			
		fileout2<<(i+1)*L<<" "<<sum_prog[i]<<" "<<err_prog[i]<<endl;

	}

	fileout2.close();




/****************************************ESERCIZIO 01.1.3**************************************************************/


	//utilizzo i 10^6 numeri casuali distribuiti uniformemente in [0, 1] generati in precedenza

	int m=100;		// sottointervalli in cui divido l'intervallo [0, 1]
	int n=10000;		// numero di elementi della sequenza di M numeri casuali che considero in ciascun intervallo temporale
	int l=int(M/n);		// numero di intervalli temporali
	
	double E=n/m;		//valore atteso in ciascun sottointervallo dopo n lanci

	double chi[l];		//stime del chi quadro in diversi intervalli temporali del generatore
	int n_i[m]={};		//vettore dei contatori che indicano quanti numeri casuali cadono in ciascuno degli m sottointervalli di [0, 1]


	
	ofstream fileout3("risultati_01.1.3.dat");


	for(int i=0; i<l; i++){

		for(int j=0; j<n; j++){

			int pos = j+i*n;
			
			for(int k=0; k<m; k++){
				
				double xmin=(double (k))/m;		//definisco gli estremi dei sottointervalli
				double xmax=(double (k+1))/m;
				
				if(rand[pos]>xmin and rand[pos]<xmax){
					n_i[k]++;			//incremento il contatore per l'intervallo considerato se il numero cade al suo interno
				}

			}

		}
	
		double chi_prog=0;

		for(int k=0; k<m; k++){

			chi_prog+=pow(n_i[k]-E,2);
			n_i[k]=0;			//azzero il contatore per l'intervallo temporale successivo

		}

		chi_prog/=E;
		chi[i]=chi_prog;

		fileout3<<i+1<<" "<<chi[i]<<endl;
	
	}

	fileout3.close();

	delete []rand;

	
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


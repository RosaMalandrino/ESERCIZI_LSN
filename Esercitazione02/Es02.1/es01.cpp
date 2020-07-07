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



/****************************************ESERCIZIO 02.1**************************************************************/


	int M=1000000;
	int N=100;
	int L=int(M/N);


	double *f_u = new double [M];
	double *f_r = new double [M];


	double I_u[N];		//stima dell'integrale con campionamento uniforme
	double I2_u[N];		//stima al quadrato con campionamento uniforme

	double sum_prog_u[N];	
	double sum2_prog_u[N];	
	double err_prog_u[N];


	double I_r[N];		//stima dell'integrale con importance sampling (in cui la distribuzione di campionamento è descritta dalla retta y=2(1-x)
	double I2_r[N];		//stima al quadrato con importance sampling

	double sum_prog_r[N];	
	double sum2_prog_r[N];	
	double err_prog_r[N];	




	for(int i=0; i<M; i++){

		double x_u=rnd.Rannyu();
		double x_r=1-sqrt(1-x_u);	//variabile distribuita secondo la retta y=2(1-x), campionata con il metodo della trasformata


		f_u[i]=M_PI/2 * cos(M_PI*x_u/2);			//integranda f(x) valutata nei punti uniformi
		f_r[i]= (M_PI/2 * cos(M_PI*x_r/2)) / (2*(1-x_r));	//f(x)/r(x) valutata nei punti campionati secondo la distribuzione r(x)


	}



	for(int i=0; i<N; i++){
	
		double sum_u=0;
		double sum_r=0;

		for(int j=0; j<L; j++){

			int k = j + i*L;
			sum_u+=f_u[k];
			sum_r+=f_r[k];

		}

		I_u[i]=sum_u/L;			//media dei valori di f(x) sui punti uniformi
		I2_u[i]=pow(I_u[i], 2);

		I_r[i]=sum_r/L;			//media dei valori di f(x)/r(x) sui punti distribuiti come r(x)
		I2_r[i]=pow(I_r[i], 2);

	
	}


	ofstream fileout("risultati_02.1.dat");

	for(int i=0; i<N; i++){
	
		for(int j=0; j<i+1; j++){

			sum_prog_u[i]+=I_u[j];	
			sum2_prog_u[i]+=I2_u[j];

			sum_prog_r[i]+=I_r[j];	
			sum2_prog_r[i]+=I2_r[j];	
		}

		sum_prog_u[i]/=(i+1);
		sum2_prog_u[i]/=(i+1);

		err_prog_u[i] = error(sum_prog_u,sum2_prog_u,i);


		sum_prog_r[i]/=(i+1);
		sum2_prog_r[i]/=(i+1);

		err_prog_r[i] = error(sum_prog_r,sum2_prog_r,i);

	
		fileout<<(i+1)*L<<" "<<sum_prog_u[i]<<" "<<err_prog_u[i]<<" "<<sum_prog_r[i]<<" "<<err_prog_r[i]<<endl;

	}

	fileout.close();


	delete []f_u;
	delete []f_r;

	
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

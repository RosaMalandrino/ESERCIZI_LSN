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



/****************************************ESERCIZIO 05.1**************************************************************/



	int M=1000000;		//numero di step
	int N=100;		//blocchi
	int L=int(M/N);		//numero di step in ciascun blocco
	

	//posizioni iniziale in unità di raggi di Bohr per campionare lo stato 100 e lo stato 210

	double x0_100=1;	// r circa 1.7 (valore atteso 1.5)
	double y0_100=1;
	double z0_100=1;

	double x0_210=3;	// r circa 5.2 (valore atteso 5)		
	double y0_210=3;
	double z0_210=3;
	

	int n_eq=200;
	double step_test;


	int dim=n_eq+M;

	double *xs=new double[dim];
	double *ys=new double[dim];
	double *zs=new double[dim];
	double *rs=new double[dim];


	for(int l=0; l<4; l++){

	
		int n_acc=0;		//passi accettati
		int n_tot=0;		//passi totali	

		ofstream eq, distr, fileout;

		string suffix;
		int orb_type;	//tipo di orbitale, =0 per 100, =1 per 210
		int step_type;	//tipo di passo, =0 se uniforme, =1 se gaussiano

		switch(l){

			case(0):			//orbitale 100, passo uniforme
				suffix="_100.dat";
				orb_type=0;
				step_type=0;
				break;
			case(1):			//orbitale 210, passo uniforme
				suffix="_210.dat";
				orb_type=1;
				step_type=0;
				break;
			case(2):			//orbitale 100, passo gaussiano
				suffix="_100_gauss.dat";
				orb_type=0;
				step_type=1;
				break;
			case(3):			//orbitale 210, passo gaussiano
				suffix="_210_gauss.dat";
				orb_type=1;
				step_type=1;
				break;
		}
		


		eq.open("equilibrazione"+suffix);
		distr.open("distribuzione"+suffix);

		if(orb_type==0){

			xs[0]=x0_100;
			ys[0]=y0_100;
			zs[0]=z0_100;
			step_test=1.1;		//step che dà un'accettazione di ~50%

		}

		if(orb_type==1){

			xs[0]=x0_210;
			ys[0]=y0_210;
			zs[0]=z0_210;
			step_test=2.8;		//step che dà un'accettazione di ~50%

		}

		if(step_type==1) step_test*=(2./3.);	//modifica dello step nel caso gaussiano per avere un'accettazione di ~50%




		for(int i=1; i<dim; i++){

			double xold=xs[i-1];
			double yold=ys[i-1];
			double zold=zs[i-1];

			
			double xtest, ytest, ztest;

			if(step_type==0){	//step uniforme

				xtest = xold + step_test*rnd.Rannyu(-1, 1);
				ytest = yold + step_test*rnd.Rannyu(-1, 1);
				ztest = zold + step_test*rnd.Rannyu(-1, 1);

			}

			if(step_type==1){	//step gaussiano

				//step_test/=2;

				xtest = xold + rnd.Gauss(0, step_test);
				ytest = yold + rnd.Gauss(0, step_test);
				ztest = zold + rnd.Gauss(0, step_test);

			}

			

			double rold = sqrt(pow(xold,2) + pow(yold,2) + pow(zold,2));
			double rtest = sqrt(pow(xtest,2) + pow(ytest,2) + pow(ztest,2));

			double cos_old=zold/rold;
			double cos_test=ztest/rtest;

			double A, ptest, pold;

			if(orb_type==0){
				ptest = exp(-2*rtest);		//trascuro i fattori di normalizzazione perché devo fare il rapporto delle due probabilità
				pold = exp(-2*rold);
			}

			if(orb_type==1){
				ptest = exp(-rtest)*pow(rtest*cos_test,2);	//trascuro i fattori di normalizzazione
				pold = exp(-rold)*pow(rold*cos_old,2);
			}

			A = ptest/pold;

			if(A>1){		//A=1, accetto sempre
				xs[i]=xtest;
				ys[i]=ytest;
				zs[i]=ztest;
				n_acc++;
			}	
			else{
				double rand=rnd.Rannyu();

				if(rand<=A){		//accetto
					xs[i]=xtest;
					ys[i]=ytest;
					zs[i]=ztest;
					n_acc++;
				}
				else{			//rigetto
					xs[i]=xold;
					ys[i]=yold;
					zs[i]=zold;
				}

			}

			n_tot++;


			rs[i]=sqrt(pow(xs[i],2) + pow(ys[i],2) + pow(zs[i],2));

			if(i<n_eq){
				eq<<rs[i]<<endl;
			}

			else	distr<<xs[i]<<" "<<ys[i]<<" "<<zs[i]<<endl;

		}

		cout<<"Rate di accettazione: "<<((double) n_acc)/n_tot<<endl;

		eq.close();
		distr.close();



		double r_ave[N];		//stima del valor medio <r>
		double r2_ave[N];		//stima del valor medio <r>^2

		double sum_prog[N]={};	
		double sum2_prog[N]={};	
		double err_prog[N];


		fileout.open("risultati_05.1" + suffix);


		for(int i=0; i<N; i++){
		
			double sum=0;

			for(int j=0; j<L; j++){

				int k = j + i*L;
				sum+=rs[n_eq+k];	//ignoro i primi n_eq punti
				
			}

			
			r_ave[i]=sum/L;			//media dei valori di r distribuiti secondo psi
			r2_ave[i]=pow(r_ave[i], 2);

		
			for(int j=0; j<i+1; j++){

				sum_prog[i]+=r_ave[j];	
				sum2_prog[i]+=r2_ave[j];
					
			}


			sum_prog[i]/=(i+1);
			sum2_prog[i]/=(i+1);

			

			err_prog[i] = error(sum_prog,sum2_prog,i);


			fileout<<(i+1)*L<<" "<<sum_prog[i]<<" "<<err_prog[i]<<endl;

		}

		fileout.close();


	}

	delete []xs;
	delete []ys;
	delete []zs;
	delete []rs;


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

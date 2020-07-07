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



/****************************************ESERCIZIO 02.2**************************************************************/


	int M=100000;		//numero di simulazioni di random walk
	int N=100;		//blocchi
	int L=int(M/N);		//numero di sottosimulazioni in ciascun blocco

	int n_steps=100;	//numero di passi del random walk


	double a=1.;

	double x_d[M]={};		//coordinate x dei punti durante le M simulazioni (caso discreto)
	double y_d[M]={};		//coordinate y dei punti durante le M simulazioni (caso discreto)
	double z_d[M]={};		//coordinate z dei punti durante le M simulazioni (caso discreto)

	double x_c[M]={};		//coordinate x dei punti durante le M simulazioni (caso continuo)
	double y_c[M]={};		//coordinate y dei punti durante le M simulazioni (caso continuo)
	double z_c[M]={};		//coordinate z dei punti durante le M simulazioni (caso continuo)



	double r_d[n_steps]={0};	//stime dello scarto quadratico medio per ciascuno degli n_steps passi	(caso discreto)
	double r2_d[n_steps]={0};	//stime al quadrato necessarie per il calcolo dell'errore		(caso discreto)
	double err_d[n_steps]={0};	//errori della stima di r per ciascuno degli n_steps passi		(caso discreto)

	double r_c[n_steps]={0};	// ""		""		""		""		""	(caso continuo)
	double r2_c[n_steps]={0};	// ""		""		""		""		""	(caso continuo)
	double err_c[n_steps]={0};	// ""		""		""		""		""	(caso continuo)

	//il primo elemento è inizializzato a 0 perché ogni random walk parte da 0, quindi la media sul passo 0 darà <r>=0; l'errore è posto anch'esso =0




	ofstream fileout("risultati_02.2.dat");


	for(int i=1; i<n_steps; i++){

		double media_d=0;
		double media2_d=0;

		double media_c=0;
		double media2_c=0;


		for(int j=0; j<N; j++){

			double sum_d=0;
			double sum_c=0;


			for(int k=0; k<L; k++){

				int pos = k + j*L;


				/************************ CASO DISCRETO *************************/


				int check=(int)rnd.Rannyu(0, 6);
			
				
				switch(check){			//spostamento di a avanti o indietro in una delle tre direzioni

					case(0):
						x_d[pos]+=a;
						break;
					
					case(1):
						x_d[pos]-=a;
						break;
					
					case(2):
						y_d[pos]+=a;
						break;
						
					case(3):
						y_d[pos]-=a;
						break;
						
					case(4):
						z_d[pos]+=a;
						break;
						
					case(5):
						z_d[pos]-=a;
						break;
						
				}
				
				sum_d+=pow(x_d[pos],2) + pow(y_d[pos],2) + pow(z_d[pos],2);		//accumulazione del valore di r^2



				/************************ CASO CONTINUO ****************************/
				

				double r=0;
				double theta=0;
				double phi=0;
				
				do{					//generazione di punti casuali all'interno della sfera di raggio a
					double x=rnd.Rannyu(-a, a);
					double y=rnd.Rannyu(-a, a);
					double z=rnd.Rannyu(-a, a);
					
					r=sqrt(pow(x,2)+pow(y,2)+pow(z,2));

					theta=acos(z/r);
					phi=atan2(y,x);
		
				}while(r>a);

				x_c[pos]+=a*sin(theta)*cos(phi);		//coordinate sulla superficie sferica di raggio a
				y_c[pos]+=a*sin(theta)*sin(phi);
				z_c[pos]+=a*cos(theta);

				sum_c+=pow(x_c[pos],2) + pow(y_c[pos],2) + pow(z_c[pos],2);		//accumulazione del valore di r^2


			}


			media_d+=sum_d/L;		//accumulazione delle N medie degli L elementi del blocco
			media2_d+=pow(sum_d/L,2);

			media_c+=sum_c/L;
			media2_c+=pow(sum_c/L,2);
		

		}


		r_d[i]=media_d/N;			//media delle medie sugli N blocchi
		r2_d[i]=media2_d/N;
		r_c[i]=media_c/N;
		r2_c[i]=media2_c/N;

		err_d[i]=sqrt((r2_d[i] - pow(r_d[i], 2))/(N-1));	//calcolo dell'errore sulla media della media su tutti gli N blocchi	
		err_c[i]=sqrt((r2_c[i] - pow(r_c[i], 2))/(N-1));


	}




	for(int i=0; i<n_steps; i++){

		//stampo la radice della media del modulo quadro

		if(i==0)	fileout<<i<<" "<<sqrt(r_d[i])<<" "<<0<<" "<<sqrt(r_c[i])<<" "<<0<<endl;		//la propagazione darebbe una divisione per 0

		else		fileout<<i<<" "<<sqrt(r_d[i])<<" "<<err_d[i]/(2*sqrt(r_d[i]))<<" "<<sqrt(r_c[i])<<" "<<err_c[i]/(2*sqrt(r_c[i]))<<endl;
		

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

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

using namespace std;


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



/****************************************ESERCIZIO 01.2**************************************************************/

	
	int n=10000;		//numero di somme da effettuare per ogni tipo di dado

	//vettori delle somme di 1, 2, 10 e 100 elementi per ogni tipo di dado

	double S1_std[n];
	double S2_std[n];
	double S10_std[n];
	double S100_std[n];

	double S1_exp[n];
	double S2_exp[n];
	double S10_exp[n];
	double S100_exp[n];

	double S1_lor[n];
	double S2_lor[n];
	double S10_lor[n];
	double S100_lor[n];
	

	ofstream std_dice("dado_standard.dat");
	ofstream exp_dice("dado_esponenziale.dat");
	ofstream lor_dice("dado_lorentziano.dat");


	for(int i=0; i<n; i++){

		double sum_std=0;		//somma di variabili casuali uniformi
		double sum_exp=0;		//somma di variabili casuali esponenziali
		double sum_lor=0;		//somma di variabili casuali lorentziane

		for(int j=0; j<100; j++){

			sum_std+=rnd.Rannyu();
			sum_exp+=rnd.Exp(1.);
			sum_lor+=rnd.Lorentz(0.,1.);


			switch(j){

				case(0):				//i-esimo risultato della somma di 1 elemento per ogni tipo di dado
					S1_std[i]=sum_std;
					S1_exp[i]=sum_exp;
					S1_lor[i]=sum_lor;
					break;

				case(1):				//i-esimo risultato della somma di 2 elementi per ogni tipo di dado
					S2_std[i]=sum_std/2;
					S2_exp[i]=sum_exp/2;
					S2_lor[i]=sum_lor/2;
					break;

				case(9):				//i-esimo risultato della somma di 10 elementi per ogni tipo di dado
					S10_std[i]=sum_std/10;
					S10_exp[i]=sum_exp/10;
					S10_lor[i]=sum_lor/10;
					break;

				case(99):				//i-esimo risultato della somma di 100 elementi per ogni tipo di dado
					S100_std[i]=sum_std/100;
					S100_exp[i]=sum_exp/100;
					S100_lor[i]=sum_lor/100;
					break;

			}		

		}

	std_dice<<S1_std[i]<<" "<<S2_std[i]<<" "<<S10_std[i]<<" "<<S100_std[i]<<endl;
	exp_dice<<S1_exp[i]<<" "<<S2_exp[i]<<" "<<S10_exp[i]<<" "<<S100_exp[i]<<endl;
	lor_dice<<S1_lor[i]<<" "<<S2_lor[i]<<" "<<S10_lor[i]<<" "<<S100_lor[i]<<endl;
		
	}

	
	std_dice.close();
	exp_dice.close();
	lor_dice.close();

	
	return 0; 
}


#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <vector>
#include <string>
#include <iomanip>

using namespace std;


double mean(vector<double> v);
double var(vector<double> v);
double AC(vector<double> v, int t);
double Error(vector<double> v, int L);


int main()
{

	string suffix;

	//Dati
	int N=500000;

	vector<double> U(N);
	vector<double> P(N);

	//Autocorrelazione

	int nmax;

	//Errore
	int stepL = 20;
	int Lmax=5000;
	int N_L = int(Lmax/stepL);
	vector<double> Ls(N_L);

	Ls[0]=10;

	for(int i=1; i<N_L; i++)	Ls[i]=Ls[i-1]+stepL;

	Ls.push_back(Lmax);

	N_L=Ls.size();

	vector<double> err_U(N_L);
	vector<double> err_P(N_L);


	for(int index=0; index<3; index++){

		
		switch(index){

			case(0):
				suffix="_solid.dat";
				nmax=200;
				break;
			case(1):
				suffix="_liquid.dat";
				nmax=400;
				break;
			case(2):
				suffix="_gas.dat";
				nmax=50;
				break;

		}

		ifstream readData;
		readData.open("risultati" + suffix);

		for (int i=0; i<N; ++i)	readData >> U[i] >> P[i];

		readData.close();

		/************* autocorrelazione **************/


		int dim=nmax;
		vector<int> n_steps(dim);

		vector<double> AC_U(dim);
		vector<double> AC_P(dim);

		for (int i=0; i<dim; ++i)	n_steps[i]=i+1;


		for(int i=0; i<dim; i++){

			int n=n_steps[i];
			AC_U[i]=AC(U, n);
			AC_P[i]=AC(P, n);

		}

		ofstream autocorr;

		autocorr.open("Analisi/AC"+suffix);

		for(int i=0; i<dim; i++)	autocorr<<n_steps[i]<<" "<<AC_U[i]<<" "<<AC_P[i]<<endl;

		autocorr.close();


		/******************* errore *********************/

		for(int i=0; i<Ls.size(); i++){

			int L=Ls[i];
			err_U[i]=Error(U, L);
			err_P[i]=Error(P, L);

		}

		ofstream err;

		err.open("Analisi/err"+suffix);

		for(int i=0; i<Ls.size(); i++)	err<<Ls[i]<<" "<<err_U[i]<<" "<<err_P[i]<<endl;

		err.close();
	

	}

	

	return 0;
}

double mean(vector<double> v){

	int size=v.size();
	double sum=0;

	for(int i=0; i<size; i++)	sum+=v[i];

	return sum/size;

}

double var(vector<double> v){

	int size=v.size();
	vector<double> v2(size);

	for(int i=0; i<size; i++)	v2[i]=v[i]*v[i];

	double m=mean(v);

	return mean(v2)-m*m;


}


double AC(vector<double> v, int t){		//t è espresso come numero di step
    
	int tmax= v.size()-1;
	vector<double> cross(tmax-t);
	vector<double> v_red(tmax-t);
	vector<double> v_shift(tmax-t);

	for(int i=0; i<tmax-t; i++){
		cross[i]=v[i]*v[i+t];
		v_red[i]=v[i];
		v_shift[i]=v[t+i];
	}

	return (mean(cross)-mean(v_red)*mean(v_shift))/var(v);

}



double Error(vector<double> v, int L){
    
	int M=v.size();
	int N=int(M/L);

	vector<double> ave(N);
	vector<double> av2(N);


	for(int i=0; i<N; i++){

		double sum = 0;

		for(int j=0; j<L; j++){
			int k = j+i*L;
			sum += v[k];
		}

		ave[i] = sum/L; 
		av2[i] = pow(ave[i],2); 

	}

	double m=mean(ave);

	return sqrt((mean(av2) - m*m)/(N-1));

}


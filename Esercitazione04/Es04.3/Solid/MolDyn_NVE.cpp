/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <iomanip>
#include "MolDyn_NVE.h"

using namespace std;

int main(){ 


	

	Input();             //Inizialization
	//int nconf = 1;


	for(int istep=1; istep <= nstep; ++istep){

		Move();           //Move particles with Verlet algorithm

		if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;

		if(istep%10 == 0){
			Measure();     			//Properties measurement
			BlockAveraging(istep/10);	//Medie a blocchi

			//ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
			//nconf += 1;
	
		}

		
	
	}



	ConfFinal();		//Write final configuration to restart
	OldFinal();		//Scrivo la penultima configurazione per il restart






	return 0;

}


void Input(void){ //Prepare all stuff for the simulation

	ifstream ReadInput,ReadConf, ReadOld;
	double ep, ek, pr, et, vir;

	cout << "Classic Lennard-Jones fluid        " << endl;
	cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
	cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
	cout << "The program uses Lennard-Jones units " << endl;

	seed = 7;    //Set seed for random numbers
	srand(seed); //Initialize random number generator
  
	ReadInput.open("input.dat"); //Read input

	ReadInput >> temp;

	ReadInput >> npart;
	cout << "Number of particles = " << npart << endl;

	ReadInput >> rho;
	cout << "Density of particles = " << rho << endl;
	vol = (double)npart/rho;
	cout << "Volume of the simulation box = " << vol << endl;
	box = pow(vol,1.0/3.0);
	cout << "Edge of the simulation box = " << box << endl;

	ReadInput >> rcut;
	ReadInput >> delta;
	ReadInput >> nstep;
	ReadInput >> iprint;
	ReadInput >> restart;	//leggo se è già stata fatta una simulazione

	cout << "The program integrates Newton equations with the Verlet method " << endl;
	cout << "Time step = " << delta << endl;
	cout << "Number of steps = " << nstep << endl;
	cout << "Restart =" << restart <<endl << endl;
	ReadInput.close();

//Prepare array for measurements
	iv = 0; //Potential energy
	ik = 1; //Kinetic energy
	ie = 2; //Total energy
	it = 3; //Temperature
	n_props = 4; //Number of observables

//Read initial configuration

	cout << "Read initial configuration from file config.0 " << endl << endl;
	ReadConf.open("config.0");

	for (int i=0; i<npart; ++i){
		ReadConf >> x[i] >> y[i] >> z[i];
		x[i] = x[i] * box;
		y[i] = y[i] * box;
		z[i] = z[i] * box;
	}

	ReadConf.close();


	if(restart==true){	//se ho fatto una simulazione in precedenza leggo anche la penultima configurazione

		cout << "Read old configuration from file old.0 " << endl << endl;
		ReadOld.open("old.0");

		for (int i=0; i<npart; ++i){
			ReadOld >> xold[i] >> yold[i] >> zold[i];
			xold[i] = xold[i] * box;
			yold[i] = yold[i] * box;
			zold[i] = zold[i] * box;
		}

		ReadOld.close();

		//x=x(t) e xold=x(t-dt)

		Move();

		//x=x(t+dt) e xold=x(t)


		double sumv2 = 0.0, fs, T;

		for(int i=0; i<npart; ++i){ //Calcolo v(t+dt/2)

			vx[i] = Pbc(x[i] - xold[i])/(delta);
			vy[i] = Pbc(y[i] - yold[i])/(delta);
			vz[i] = Pbc(z[i] - zold[i])/( delta);

			sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];

		}

		sumv2 /= (double)npart;

		T=sumv2/3;

		fs=sqrt(temp/T);	//trovo il fattore di scala affinché le velocità rispecchino la temperatura T da simulare


		for (int i=0; i<npart; ++i){

			vx[i] *= fs;
			vy[i] *= fs;
			vz[i] *= fs;

			xold[i] = Pbc(x[i] - vx[i] * delta);
			yold[i] = Pbc(y[i] - vy[i] * delta);
			zold[i] = Pbc(z[i] - vz[i] * delta);
		
		}

		

	//Alla fine x = r(t+deltat) e xold=rnew


	}


	if(restart==false){
	//Prepare initial velocities
		cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
		double sumv[3] = {0.0, 0.0, 0.0};

		for (int i=0; i<npart; ++i){
			vx[i] = rand()/double(RAND_MAX) - 0.5;
			vy[i] = rand()/double(RAND_MAX) - 0.5;
			vz[i] = rand()/double(RAND_MAX) - 0.5;

			sumv[0] += vx[i];
			sumv[1] += vy[i];
			sumv[2] += vz[i];
		}

		for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;

		double sumv2 = 0.0, fs;

		for (int i=0; i<npart; ++i){
			vx[i] = vx[i] - sumv[0];
			vy[i] = vy[i] - sumv[1];
			vz[i] = vz[i] - sumv[2];

			sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
		}
		sumv2 /= (double)npart;

		fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 

		for (int i=0; i<npart; ++i){

			vx[i] *= fs;
			vy[i] *= fs;
			vz[i] *= fs;

			xold[i] = Pbc(x[i] - vx[i] * delta);
			yold[i] = Pbc(y[i] - vy[i] * delta);
			zold[i] = Pbc(z[i] - vz[i] * delta);
		
		}
	}

   
	return;

}


void Move(void){ //Move particles with Verlet algorithm

	double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

	for(int i=0; i<npart; ++i){ //Force acting on particle i
	
		fx[i] = Force(i,0);
		fy[i] = Force(i,1);
		fz[i] = Force(i,2);

	}

	for(int i=0; i<npart; ++i){ //Verlet integration scheme

		xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
		ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
		znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

		vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
		vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
		vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

		xold[i] = x[i];
		yold[i] = y[i];
		zold[i] = z[i];

		x[i] = xnew;
		y[i] = ynew;
		z[i] = znew;

	}

	return;

}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  
	double f=0.0;
	double dvec[3], dr;

	for (int i=0; i<npart; ++i){

		if(i != ip){

			dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
			dvec[1] = Pbc( y[ip] - y[i] );
			dvec[2] = Pbc( z[ip] - z[i] );

			dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
			dr = sqrt(dr);

			if(dr < rcut){
				f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
			}
		
		}
  
	}
  
	return f;

}

void Measure(){ //Properties measurement
  
	int bin;
	double v, t, vij;
	double dx, dy, dz, dr;
	ofstream Epot, Ekin, Etot, Temp;

	Epot.open("output_epot.dat",ios::app);
	Ekin.open("output_ekin.dat",ios::app);
	Temp.open("output_temp.dat",ios::app);
	Etot.open("output_etot.dat",ios::app);

	v = 0.0; //reset observables
	t = 0.0;

//cycle over pairs of particles
	for (int i=0; i<npart-1; ++i){

		for (int j=i+1; j<npart; ++j){

			dx = Pbc( xold[i] - xold[j] );	//correzione spero giusta x-----> xold, perché dopo aver fatto Move xold=x(t), mentre x=x(t+deltat),
			dy = Pbc( yold[i] - yold[j] );	//mentre Ekin è al tempo t
			dz = Pbc( zold[i] - zold[j] );

			dr = dx*dx + dy*dy + dz*dz;
			dr = sqrt(dr);

			if(dr < rcut){
				vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);

				//Potential energy
				v += vij;
			}
		
		}          
 
	}

//Kinetic energy
	for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
	stima_pot = v/(double)npart; //Potential energy per particle
	stima_kin = t/(double)npart; //Kinetic energy per particle
	stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
	stima_etot = (t+v)/(double)npart; //Total energy per particle

	Epot << stima_pot  << endl;
	Ekin << stima_kin  << endl;
	Temp << stima_temp << endl;
	Etot << stima_etot << endl;

	Epot.close();
	Ekin.close();
	Temp.close();
	Etot.close();

	return;

}



void BlockAveraging(int istep){ //Medie a blocchi, passo il punto della simulazione in cui sono, diviso per ogni quante volte misuro


	int L=int(nstep/(10*nblocks));	//numero di passi per blocco (considerando che faccio una misura ogni 10 step)
  
	if(istep%(L)==1){	//accumulatori =0 all'inizio del blocco

		sum_pot=0;
		sum_kin=0;
		sum_etot=0;
		sum_temp=0;

	}


	sum_pot+=stima_pot;
	sum_kin+=stima_kin;
	sum_etot+=stima_etot;
	sum_temp+=stima_temp;



	

	if(istep%(L)==0){	//è finito il blocco posso fare la media

		int i=istep/(L)-1;		//indice del vettore delle medie in ciascun blocco, il rapporto va da 1 a nblocks, sottraggo 1

		ave_pot[i] = sum_pot/L;			
		ave2_pot[i] = pow(ave_pot[i],2);

		ave_kin[i] = sum_kin/L;			
		ave2_kin[i] = pow(ave_kin[i],2);

		ave_etot[i] = sum_etot/L;			
		ave2_etot[i] = pow(ave_etot[i],2);	
		
		ave_temp[i] = sum_temp/L;			
		ave2_temp[i] = pow(ave_temp[i],2);

		
		sum_prog_pot[i]=0;
		sum_prog_kin[i]=0;
		sum_prog_etot[i]=0;
		sum_prog_temp[i]=0;		

		
		for(int j=0; j<i+1; j++){

			sum_prog_pot[i]+=ave_pot[j];	
			sum2_prog_pot[i]+=ave2_pot[j];

			sum_prog_kin[i]+=ave_kin[j];	
			sum2_prog_kin[i]+=ave2_kin[j];

			sum_prog_etot[i]+=ave_etot[j];	
			sum2_prog_etot[i]+=ave2_etot[j];

			sum_prog_temp[i]+=ave_temp[j];	
			sum2_prog_temp[i]+=ave2_temp[j];

		}


		sum_prog_pot[i]/=(i+1);		
		sum2_prog_pot[i]/=(i+1);

		sum_prog_kin[i]/=(i+1);		
		sum2_prog_kin[i]/=(i+1);

		sum_prog_etot[i]/=(i+1);		
		sum2_prog_etot[i]/=(i+1);	

		sum_prog_temp[i]/=(i+1);		
		sum2_prog_temp[i]/=(i+1);


		const int wd=15;
		const int cfr=9;	//cifre da stampare

		ofstream Epot_ave, Ekin_ave, Etot_ave, Temp_ave;

		Epot_ave.open("ave_epot.dat",ios::app);
		Ekin_ave.open("ave_ekin.dat",ios::app);
		Etot_ave.open("ave_etot.dat",ios::app);
		Temp_ave.open("ave_temp.dat",ios::app);

		int step_blk=nstep/nblocks;


		Epot_ave << setw(wd) << (i+1)*step_blk << setw(wd) << setprecision(cfr) << sum_prog_pot[i] << setw(wd) << setprecision(cfr) << error(sum_prog_pot,sum2_prog_pot,i)<<endl;
		Ekin_ave << setw(wd) << (i+1)*step_blk << setw(wd) << setprecision(cfr) << sum_prog_kin[i] << setw(wd) << setprecision(cfr) << error(sum_prog_kin,sum2_prog_kin,i)<<endl;
		Etot_ave << setw(wd) << (i+1)*step_blk << setw(wd) << setprecision(cfr) << sum_prog_etot[i]<< setw(wd) << setprecision(cfr) << error(sum_prog_etot,sum2_prog_etot,i)<<endl;
		Temp_ave << setw(wd) << (i+1)*step_blk << setw(wd) << setprecision(cfr) << sum_prog_temp[i]<< setw(wd) << setprecision(cfr) << error(sum_prog_temp,sum2_prog_temp,i)<<endl;
		

		Epot_ave.close();
		Ekin_ave.close();
		Temp_ave.close();
		Etot_ave.close();


	}


	return;

}




void ConfFinal(void){ //Write final configuration
  
	ofstream WriteConf;

	cout << "Print final configuration to file config.final " << endl << endl;
	WriteConf.open("config.final");

	for (int i=0; i<npart; ++i){

		WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  
	}

	WriteConf.close();

	return;
}


void OldFinal(void){ //Write final old configuration
  
	ofstream WriteOld;

	cout << "Print final old configuration to file old.final " << endl << endl;
	WriteOld.open("old.final");

	for (int i=0; i<npart; ++i){

		WriteOld << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
  
	}

	WriteOld.close();

	return;
}


void ConfXYZ(int nconf){ //Write configuration in .xyz format
  
	ofstream WriteXYZ;

	WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
	WriteXYZ << npart << endl;
	WriteXYZ << "This is only a comment!" << endl;
  
	for (int i=0; i<npart; ++i){

		WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  
	}

	WriteXYZ.close();

}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    
	return r - box * rint(r/box);

}




double error(double av[], double av2[], int n){


	if (n == 0){
		return 0;
	}
	else{
		 return sqrt((av2[n]- (av[n])*(av[n]))/n);			
	}		


}



/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

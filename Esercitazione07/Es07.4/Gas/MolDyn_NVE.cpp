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
	int nconf = 1;


	for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
	{
		Reset(iblk);   //Reset block averages
		for(int istep=1; istep <= nstep; ++istep)
		{

			int k=(iblk-1)*nstep+istep;

			Move();

			if(k%iprint == 0) cout << "Number of time-steps: " << k << endl;

			
			if(k%10 == 0){
				Measure();
				Accumulate(); //Update block averages
				//ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
				nconf += 1;
			}
		}
		Averages(iblk);   //Print results for current block
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
	ReadInput >> nblk;
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

//measurement of g(r)
	igofr = n_props;
	nbins = 100;
	n_props = n_props + nbins;
	bin_size = (box/2.0)/(double)nbins;

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

		//update of the histogram of g(r)

			if(dr < box/2.0)
			{
			bin = int(dr/bin_size);
			walker[igofr+bin]+=2;

			}
			


			if(dr < rcut){
				vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);

				//Potential energy
				v += vij;
			}
		
		}          
 
	}

//Kinetic energy
	for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
	walker[iv] = v/(double)npart; //Potential energy per particle
	walker[ik] = t/(double)npart; //Kinetic energy per particle
	walker[it] = (2.0 / 3.0) * t/(double)npart; //Temperature
	walker[ie] = (t+v)/(double)npart; //Total energy per particle

	Epot << walker[iv]  << endl;
	Ekin << walker[ik]  << endl;
	Temp << walker[it] << endl;
	Etot << walker[ie] << endl;


	Epot.close();
	Ekin.close();
	Temp.close();
	Etot.close();


	for (int k=igofr; k<igofr+nbins; ++k){

		int j = k-igofr;
		double dV=(4*pi/3)*(pow((j+1)*bin_size, 3) - pow(j*bin_size, 3));
		//dV=1.0;
		walker[k]/=rho*npart*dV;

		//cout<<walker[k]<<endl;

	}

	return;

}


void Reset(int iblk) //Reset block averages
{
   
	if(iblk == 1)
	{
		for(int i=0; i<n_props; ++i)
		{
			glob_av[i] = 0;
			glob_av2[i] = 0;
		}
	}

	for(int i=0; i<n_props; ++i)
	{
		blk_av[i] = 0;
	}
	blk_norm = 0;
	attempted = 0;
	accepted = 0;
}


void Accumulate(void) //Update block averages
{

	for(int i=0; i<n_props; ++i)
	{
		blk_av[i] = blk_av[i] + walker[i];
	}
	blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
    
	double r, gdir;
	
	const int wd=12;
    
	cout << "Block number " << iblk << endl;
	//cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
	ofstream Epot_ave, Ekin_ave, Etot_ave, Temp_ave;
	ofstream Gofr, Gave;

	Epot_ave.open("ave_epot.dat",ios::app);
	Ekin_ave.open("ave_ekin.dat",ios::app);
	Etot_ave.open("ave_etot.dat",ios::app);
	Temp_ave.open("ave_temp.dat",ios::app);

	Gofr.open("output.gofr.0",ios::app);
	Gave.open("output.gave.0",ios::app);
    
	//stima_pot = blk_av[iv]/blk_norm/(double)npart + vtail; //Potential energy
	stima_pot = blk_av[iv]/blk_norm/(double)npart; //Potential energy
	glob_av[iv] += stima_pot;
	glob_av2[iv] += stima_pot*stima_pot;
	err_pot=Error(glob_av[iv],glob_av2[iv],iblk);
    
	stima_kin = blk_av[ik]/blk_norm/(double)npart; //Kinetic energy
	glob_av[ik] += stima_kin;
	glob_av2[ik] += stima_kin*stima_kin;
	err_kin=Error(glob_av[ik],glob_av2[ik],iblk);

	stima_etot = blk_av[ie]/blk_norm/(double)npart; //Total energy
	glob_av[ie] += stima_etot;
	glob_av2[ie] += stima_etot*stima_etot;
	err_etot=Error(glob_av[ie],glob_av2[ie],iblk);

	stima_temp = blk_av[it]/blk_norm/(double)npart; //Temperature
	glob_av[it] += stima_temp;
	glob_av2[it] += stima_temp*stima_temp;
	err_temp=Error(glob_av[it],glob_av2[it],iblk);


	for (int k=igofr; k<igofr+nbins; ++k){
		
		//int j = k-igofr;
		stima_gdir = blk_av[k]/blk_norm;
		glob_av[k] += stima_gdir;
		glob_av2[k] += stima_gdir*stima_gdir;
		//err_gdir=Error(glob_av[k],glob_av2[k],iblk);

	}

//Potential energy per particle
	Epot_ave << setw(wd) << iblk <<  setw(wd) << stima_pot << setw(wd) << glob_av[iv]/(double)iblk << setw(wd) << err_pot << endl;
//Kinetic energy per particle
	Ekin_ave << setw(wd) << iblk <<  setw(wd) << stima_kin << setw(wd) << glob_av[ik]/(double)iblk << setw(wd) << err_kin << endl;
//Total energy per particle
	Etot_ave << setw(wd) << iblk <<  setw(wd) << stima_etot<< setw(wd) << glob_av[ie]/(double)iblk << setw(wd) << err_etot<< endl;
//Pressure
	Temp_ave << setw(wd) << iblk <<  setw(wd) << stima_temp<< setw(wd) << glob_av[it]/(double)iblk << setw(wd) << err_temp<< endl;


//g(r)

	Gofr<< iblk <<  setw(wd);
	for(int k=igofr; k<igofr+nbins; ++k)	Gofr << glob_av[k]/(double)iblk << setw(wd);	//Gofr << glob_av[k]/(double)iblk << "\t";
	Gofr << endl;

	if(iblk == nblk-1){
		
		for(int k=igofr; k<igofr+nbins; ++k){

			err_gdir=Error(glob_av[k],glob_av2[k],iblk);
			Gave << (k-igofr + 0.5)*bin_size <<  setw(wd) << glob_av[k]/(double)iblk << setw(wd) << err_gdir << endl;

		}
	}

	cout << "----------------------------" << endl << endl;

	Epot_ave.close();
	Ekin_ave.close();
	Etot_ave.close();
	Temp_ave.close();

	Gofr.close();
	Gave.close();
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




double Error(double sum, double sum2, int iblk)
{
	if( iblk == 1 ) return 0.0;
	else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
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

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
//parameters, observables
const int m_props=4;
int n_props;
int iv,ik,it,ie;
double stima_pot, stima_kin, stima_etot, stima_temp;

// averages
double acc,att;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep, iprint, seed;
double delta;
bool restart; //vedo se ricomincio o no

//block averaging
double sum_pot, sum_kin, sum_etot, sum_temp;
const int nblocks=100;
double ave_pot[nblocks], ave2_pot[nblocks], sum_prog_pot[nblocks], sum2_prog_pot[nblocks];
double ave_kin[nblocks], ave2_kin[nblocks], sum_prog_kin[nblocks], sum2_prog_kin[nblocks];
double ave_etot[nblocks], ave2_etot[nblocks], sum_prog_etot[nblocks], sum2_prog_etot[nblocks];			
double ave_temp[nblocks], ave2_temp[nblocks], sum_prog_temp[nblocks], sum2_prog_temp[nblocks];


//functions
void Input(void);
void Move(void);
void ConfFinal(void);
void OldFinal(void);
void ConfXYZ(int);
void Measure(void);
void BlockAveraging(int);
double Force(int, int);
double Pbc(double);

double error(double [], double [], int);
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

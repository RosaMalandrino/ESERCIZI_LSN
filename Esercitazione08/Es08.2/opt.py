import numpy as np
import os

np.random.seed(3)

mu0 = 2.0
sigma0 = 1.0

step_mu = 0.1
step_sigma= 0.1

n_h=10		#numero di suddivisioni delle temperature alte
n_l=10		#numero di suddivisioni delle temperature basse

n=n_h+n_l	#numero di suddivisioni della temperatura

N=50		#numero di step per ciascuna temperatura

mus= np.zeros(N*n+1)
sigmas = np.zeros(N*n+1)
Es = np.zeros(N*n+1)


exists = os.path.exists('ottimizzazione.dat')

if exists:
	os.system('rm ottimizzazione.dat')

par1=str(mu0)
par2=str(sigma0)

os.system('./VarMonteCarlo.exe '+ par1 + ' ' + par2)

E0 = np.loadtxt("ottimizzazione.dat", usecols=(0), delimiter=' ', unpack=True)


mus[0]=mu0
sigmas[0] = sigma0
Es[0]=E0



T_l=0.1
Tmax=1 + T_l


Ts_h=np.array([Tmax - i*(Tmax-T_l)/n_h for i in range(n_h)])	#n_h temperature tra 1.1 e 0.1
Ts_l=np.array([T_l - i*T_l/n_l for i in range(n_l)])		#n_l temperature tra 0.1 e 0
Ts=np.concatenate((Ts_h, Ts_l))


#print(Ts)

for j in range(n):

	T=Ts[j]
	beta=1/T

	n_acc=0
	n_tot=0


	if j==int(0.5*n):		#a metà dei campioni di T faccio step più piccoli
		step_mu = 0.05
		step_sigma= 0.05


	for k in range(N):

		i=j*N + k +1

		mu_old=mus[i-1]
		sigma_old=sigmas[i-1]
		E_old=Es[i-1]

		mu_test = mu_old + step_mu*np.random.uniform(-1, 1)
		sigma_test = sigma_old + step_sigma*np.random.uniform(-1, 1)

		if mu_test < 0:
			mu_test = -mu_test

		if sigma_test<0:
			sigma_test = -sigma_test

		par1=str(mu_test)
		par2=str(sigma_test)

		os.system('./VarMonteCarlo.exe '+ par1 + ' ' + par2)

		E_test = np.loadtxt("ottimizzazione.dat", usecols=(0), delimiter=' ', unpack=True, skiprows=i)	


		################################################# METROPOLIS ############################################################


		p=np.exp(-beta*(E_test - E_old))

		r=np.random.uniform()

		if r<=p:	#accetto
			mus[i]=mu_test
			sigmas[i]=sigma_test
			Es[i]=E_test
			n_acc = n_acc + 1

		else:
			mus[i]=mu_old
			sigmas[i]=sigma_old
			Es[i]=E_old

		n_tot = n_tot + 1


	print("\nTemperatura: ", round(T, 3))
	print("Rate di accettazione dell'ottimizzazione: ", n_acc/n_tot, "\n")
		

	

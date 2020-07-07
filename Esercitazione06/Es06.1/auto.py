import numpy as np
import os

Tmax=2.0
Tmin=0.5
n=21

dT=(Tmax-Tmin)/(n-1)

campo=0.0
metro=1

os.system("./clean.sh")
os.system("rm -rf risultati_*")
os.system("rm -rf equilibrazione.dat")



for j in range(2):

	if j==0:
		metro=1

	else:	metro=0

	restart=0

	for i in range(2*n):

		if i>1:
			restart=1

		if i%2==0:
			T = Tmax - (i/2)*dT	#sovrascivo la temperatura per studiare l'andamento	
			campo = 0.0

		elif i%2==1:
			campo = 0.02


		file = open("input.dat", "r")
		lines = file.readlines()
		lines[0] = str(T)+"\n"
		lines[3] = str(campo)+"\n"
		lines[4] = str(metro)+"\n"
		lines[7] = str(restart)+"\n"

		file = open("input.dat", "w")
		file.writelines(lines)
		file.close()


		os.system("./Monte_Carlo_ISING_1D.exe")

		
		if i%2==0:
			os.system("cp config.final config.restart")

		os.system("./clean.sh")


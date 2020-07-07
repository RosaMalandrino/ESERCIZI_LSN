import os

Tmin=0.2
Tmax=4
n=20
dT=(Tmax-Tmin)/(n-1)

Ts = [Tmin + i*dT for i in range(n)]

os.system("cp input.pimc input.dat")
os.system("rm prob_pimc_*")


for i in range(n):

	T=round(Ts[i],1)

	file = open("input.dat", "r")
	lines = file.readlines()
	lines[1] = "temperature\t\t\t\t"+str(T)+"\n"
	

	file = open("input.dat", "w")
	file.writelines(lines)
	file.close()

	print("Temperature: ", T, "\n")
	os.system("./qmc1d")

	print("\n")

	os.system("cp probability.dat prob_pimc_"+str(T)+".dat")


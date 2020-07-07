CC = g++
CFLAGS = -Wall -O3 --std=c++11

VarMonteCarlo.exe : VarMonteCarlo.o random.o
	$(CC) random.o VarMonteCarlo.o -o VarMonteCarlo.exe
VarMonteCarlo.o : VarMonteCarlo.cpp
	$(CC) -c VarMonteCarlo.cpp -o VarMonteCarlo.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o VarMonteCarlo.exe seed.out

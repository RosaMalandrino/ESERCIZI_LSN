CC = g++
CFLAGS = -Wall -O3 --std=c++11

es01.exe : es01.o random.o genetic.o
	$(CC) random.o genetic.o es01.o -o es01.exe
es01.o : es01.cpp
	$(CC) -c es01.cpp -o es01.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
genetic.o : genetic.cpp genetic.h random.h
	$(CC) -c genetic.cpp -o genetic.o $(CFLAGS)

clean :
	rm *.o es01.exe seed.out

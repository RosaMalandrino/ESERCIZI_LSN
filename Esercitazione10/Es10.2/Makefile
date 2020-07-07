CC = mpicxx
CFLAGS = -Wall -O3 --std=c++11

es02.exe : es02.o random.o genetic.o
	$(CC) random.o genetic.o es02.o -o es02.exe
es02.o : es02.cpp
	$(CC) -c es02.cpp -o es02.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
genetic.o : genetic.cpp genetic.h random.h
	$(CC) -c genetic.cpp -o genetic.o $(CFLAGS)

clean :
	rm *.o es02.exe seed.out

CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : Esercizio9_1.o lib.o random.o TSP_GA.o
	$(CC) lib.o random.o TSP_GA.o Esercizio9_1.o -o Esercizio9_1.exe
main.o : Esercizio9_1.cpp
	$(CC) -c Esercizio9_1.cpp -o Esercizio9_1.o $(CFLAGS)
random.o : ../../Librerie/random.cpp ../../Librerie/random.h
	$(CC) -c ../../Librerie/random.cpp -o random.o $(CFLAGS)
lib.o : ../../Librerie/lib.C ../../Librerie/lib.h
	$(CC) -c ../../Librerie/lib.C -o lib.o $(CFLAGS)
TSP_GA.o : TSP_GA.cpp TSP_GA.h
	$(CC) -c TSP_GA.cpp -o TSP_GA.o $(CFLAGS)
clean :
	rm *.o *.txt Esercizio9_1.exe seed.out
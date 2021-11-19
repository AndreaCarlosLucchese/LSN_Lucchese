CC = g++
CFLAGS = -Wall -O3 --std=c++11
CP= cp
main.exe : Esercizio5_1.o lib.o random.o hydro.o
	$(CC) lib.o random.o hydro.o Esercizio5_1.o -o Esercizio5_1.exe
main.o : Esercizio5_1.cpp
	$(CC) -c Esercizio5_1.cpp -o Esercizio5_1.o $(CFLAGS)
random.o : ../../Librerie/random.cpp ../../Librerie/random.h
	$(CC) -c ../../Librerie/random.cpp -o random.o $(CFLAGS)
lib.o : ../../Librerie/lib.C ../../Librerie/lib.h
	$(CC) -c ../../Librerie/lib.C -o lib.o $(CFLAGS)
hydro.o : hydro.c hydro.h
	$(CC) -c hydro.c -o hydro.o $(CFLAGS)
clean:
	rm *.o *.txt *.dat Esercizio5_1.exe seed.out
	make
Test: Equilibration_100.dat
	$(CP) Equilibration_100.dat Test/Equilibration_100.dat
	$(CP) Equilibration_210.dat Test/Equilibration_210.dat
CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : Esercizio8_1.o lib.o random.o Quantum.o
	$(CC) lib.o random.o Quantum.o Esercizio8_1.o -o Esercizio8_1.exe
main.o : Esercizio8_1.cpp
	$(CC) -c Esercizio8_1.cpp -o Esercizio8_1.o $(CFLAGS)
random.o : ../../Librerie/random.cpp ../../Librerie/random.h
	$(CC) -c ../../Librerie/random.cpp -o random.o $(CFLAGS)
lib.o : ../../Librerie/lib.C ../../Librerie/lib.h
	$(CC) -c ../../Librerie/lib.C -o lib.o $(CFLAGS)
Quantum.o : Quantum.cpp Quantum.h
	$(CC) -c Quantum.cpp -o Quantum.o $(CFLAGS)
clean :
	rm *.o *.txt main.exe seed.out

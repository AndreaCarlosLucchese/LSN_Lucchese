CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : Esercizio3_1.o lib.o random.o
	$(CC) lib.o random.o Esercizio3_1.o -o Esercizio3_1.exe
main.o : Esercizio3_1.cpp
	$(CC) -c Esercizio3_1.cpp -o Esercizio3_1.o $(CFLAGS)
random.o : ../../Librerie/random.cpp ../../Librerie/random.h
	$(CC) -c ../../Librerie/random.cpp -o random.o $(CFLAGS)
lib.o : ../../Librerie/lib.C ../../Librerie/lib.h
	$(CC) -c ../../Librerie/lib.C -o lib.o $(CFLAGS)
clean :
	rm *.o *.txt main.exe seed.out 

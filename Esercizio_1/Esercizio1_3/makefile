CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : Esercizio1_3.o lib.o random.o
	$(CC) lib.o random.o Esercizio1_3.o -o Esercizio1_3.exe
main.o : Esercizio1_2.cpp
	$(CC) -c Esercizio1_3.cpp -o Esercizio1_3.o $(CFLAGS)
random.o : ../../Librerie/random.cpp ../../Librerie/random.h
	$(CC) -c ../../Librerie/random.cpp -o random.o $(CFLAGS)
lib.o : ../../Librerie/lib.C ../../Librerie/lib.h
	$(CC) -c ../../Librerie/lib.C -o lib.o $(CFLAGS)
clean :
	rm *.o *.txt main.exe seed.out 

CC = mpic++
CFLAGS = -Wall -O3 --std=c++11

main.exe : Esercizio10_2.o lib.o random.o TSP_SA.o
	$(CC) lib.o random.o TSP_SA.o Esercizio10_2.o -o Esercizio10_2.exe
main.o : Esercizio10_2.cpp
	$(CC) -c Esercizio10_2.cpp -o Esercizio10_2.o $(CFLAGS)
random.o : ../../Librerie/random.cpp ../../Librerie/random.h
	$(CC) -c ../../Librerie/random.cpp -o random.o $(CFLAGS)
lib.o : ../../Librerie/lib.C ../../Librerie/lib.h
	$(CC) -c ../../Librerie/lib.C -o lib.o $(CFLAGS)
TSP_SA.o : TSP_SA.cpp TSP_SA.h
	$(CC) -c TSP_SA.cpp -o TSP_SA.o $(CFLAGS)
clean :
	rm *.o *.txt Esercizio10_1.exe seed.out
CC = g++
CFLAGS = -Wall -O3 --std=c++11
CP =cp
main.exe : MolDyn_NVE.o 
	$(CC)  MolDyn_NVE.o -o Esercizio4_1.exe
main.o : MolDyn_NVE.cpp
	$(CC) -c MolDyn_NVE.cpp -o MolDyn_NVE.o $(CFLAGS)
clean :
	rm *.out main.exe seed.out
Fcc : config.fcc config.0
	$(CP) config.fcc config.0
Conf: old.final old.0 config.final config.0
	$(CP) old.final old.0
	$(CP) config.final config.0
Solid: argon.solid input.dat #Copio configurazione argon solido
	$(CP) argon.solid input.dat	
Liquid: argon.liquid input.dat #Copio configurazione argon liquido
	$(CP) argon.liquid input.dat
Gas: argon.gas input.dat #Copio configurazione argon liquido
	$(CP) argon.gas input.dat
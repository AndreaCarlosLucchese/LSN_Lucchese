CC = g++
CFLAGS = -Wall -O3 --std=c++11
CP =cp
main.exe : MolDyn_NVE.o 
	$(CC)  MolDyn_NVE.o -o MolDyn_NVE.exe
main.o : MolDyn_NVE.cpp
	$(CC) -c MolDyn_NVE.cpp -o MolDyn_NVE.o $(CFLAGS)
clean :
	rm *.o  *.out main.exe seed.out
Fcc : config.fcc config.0 #Copio la configurazione FCC in quella di partenza
	$(CP) config.fcc config.0
Conf: old.final old.0 config.final config.0 # Copio la configurazione Finale in quella di partenza
	$(CP) old.final old.0
	$(CP) config.final config.0
Solid: input.solid input.dat #Copio i dati del solido nell'input
	$(CP) input.solid input.dat	
Liquid: input.liquid input.dat #Copio i dati del liquido nell'input
	$(CP) input.liquid input.dat

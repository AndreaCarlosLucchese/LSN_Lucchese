CC = g++
CFLAGS = -Wall -O3 --std=c++11
CP =cp
main.exe : MolDyn_NVE.o 
	$(CC)  MolDyn_NVE.o -o Esercizio4_1.exe
main.o : MolDyn_NVE.cpp
	$(CC) -c MolDyn_NVE.cpp -o MolDyn_NVE.o $(CFLAGS)
clean :
	rm *.out Esercizio4_1.exe *.o
	make
Fcc : config.fcc config.0
	$(CP) config.fcc config.0
Conf: old.final old.0 config.final config.0
	$(CP) old.final old.0
	$(CP) config.final config.0
Solid: argon.solid input.dat #Copio configurazione argon solido
	$(CP) argon.solid input.dat	
	$(CP) ave_ekin.out Solid_result/ave_ekin.out
	$(CP) ave_epot.out Solid_result/ave_epot.out
	$(CP) ave_etot.out Solid_result/ave_etot.out
	$(CP) ave_temp.out Solid_result/ave_temp.out
Liquid: argon.liquid input.dat #Copio configurazione argon liquido
	$(CP) argon.liquid input.dat
	$(CP) ave_ekin.out Liquid_result/ave_ekin.out
	$(CP) ave_epot.out Liquid_result/ave_epot.out
	$(CP) ave_etot.out Liquid_result/ave_etot.out
	$(CP) ave_temp.out Liquid_result/ave_temp.out
Gas: argon.gas input.dat #Copio configurazione argon liquido
	$(CP) argon.gas input.dat
	$(CP) ave_ekin.out Gas_result/ave_ekin.out
	$(CP) ave_epot.out Gas_result/ave_epot.out
	$(CP) ave_etot.out Gas_result/ave_etot.out
	$(CP) ave_temp.out Gas_result/ave_temp.out
	
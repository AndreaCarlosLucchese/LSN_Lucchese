CC = g++
CFLAGS = -Wall -O3 --std=c++11
CP =cp
main.exe : MolDyn_NVE.o 
	$(CC)  MolDyn_NVE.o -o MolDyn_NVE.exe
main.o : MolDyn_NVE.cpp
	$(CC) -c MolDyn_NVE.cpp -o MolDyn_NVE.o $(CFLAGS)
clean :
	rm *.o *.out MolDyn_NVE.exe 
	make
Fcc : config.fcc config.0 #Copio la configurazione FCC in quella di partenza
	$(CP) config.fcc config.0
Conf: old.final old.0 config.final config.0 # Copio la configurazione Finale in quella di partenza
	$(CP) old.final old.0
	$(CP) config.final config.0
Solid: input.solid input.dat #Copio i dati del solido nell'input
	$(CP) input.solid input.dat	
	$(CP) ave_ekin.out Solid_result/ave_ekin.out
	$(CP) ave_epot.out Solid_result/ave_epot.out
	$(CP) ave_etot.out Solid_result/ave_etot.out
	$(CP) ave_temp.out Solid_result/ave_temp.out
Liquid: input.liquid input.dat #Copio i dati del liquido nell'input
	$(CP) input.liquid input.dat
	$(CP) ave_ekin.out Liquid_result/ave_ekin.out
	$(CP) ave_epot.out Liquid_result/ave_epot.out
	$(CP) ave_etot.out Liquid_result/ave_etot.out
	$(CP) ave_temp.out Liquid_result/ave_temp.out
Equilibration: output_ekin.dat 
	$(CP) output_ekin.dat Equilibration/output_ekin.dat
	$(CP) output_epot.dat Equilibration/output_epot.dat
	$(CP) output_etot.dat Equilibration/output_etot.dat
	$(CP) output_temp.dat Equilibration/output_temp.dat
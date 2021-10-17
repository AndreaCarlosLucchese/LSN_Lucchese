CC = g++
CFLAGS = -Wall -O3 --std=c++11
CP =cp
main.exe : MolDyn_NVE.o 
	$(CC)  MolDyn_NVE.o -o MolDyn_NVE.exe
main.o : MolDyn_NVE.cpp
	$(CC) -c MolDyn_NVE.cpp -o MolDyn_NVE.o $(CFLAGS)
clean :
	rm *.o *.out MolDyn_NVE.exe
	$(CP) config.final config.0
	$(CP) old.final old.0
	make
trash :
	rm *.o *.out
esegui:
	./Monte_Carlo_NVT.exe
Fcc: config.fcc config.0
	$(CP) config.fcc config.0
Conf:  config.final config.0
	$(CP) config.final config.0
	$(CP) old.final old.0
Solid: argon.solid input.dat
	$(CP) argon.solid input.dat	
	$(CP) ave_pres.out Solid_result/ave.pres.out
	$(CP) ave_epot.out Solid_result/ave.epot.out
	$(CP) output_gave.out Solid_result/output_gave.out 
Liquid: argon.liquid input.dat
	$(CP) argon.liquid input.dat
	$(CP) ave_pres.out Liquid_result/ave.pres.out
	$(CP) ave_epot.out Liquid_result/ave.epot.out
	$(CP) output_gave.out Liquid_result/output_gave.out 
Gas: argon.gas input.dat
	$(CP) argon.gas input.dat
	$(CP) ave_pres.out Gas_result/ave.pres.out
	$(CP) ave_epot.out Gas_result/ave.epot.out
	$(CP) output_gave.out Gas_result/output_gave.out 

#identico a quello del MonteCarlo
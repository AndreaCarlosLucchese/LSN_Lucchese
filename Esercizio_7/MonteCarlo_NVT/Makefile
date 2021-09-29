CC = g++
CFLAGS = -Wall -O3 --std=c++11
CP =cp

Monte_Carlo_NVT.exe : Monte_Carlo_NVT.o random.o
	$(CC) random.o Monte_Carlo_NVT.o -o Monte_Carlo_NVT.exe
Monte_Carlo_NVT.o : Monte_Carlo_NVT.cpp
	$(CC) -c Monte_Carlo_NVT.cpp -o Monte_Carlo_NVT.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o *.0 *.txt Monte_Carlo_NVT.exe 
	$(CP) config.final config.0
	make
esegui:
	./Monte_Carlo_NVT.exe
Fcc: config.fcc config.0 #Copio configurazione FCC
	$(CP) config.fcc config.0
Conf:  config.final config.0 #Copio configurazione Finae
	$(CP) config.final config.0
Solid: input.solid input.datm # Copio i dati del solido in quelli di input e copio gli output nella relativa cartella
	$(CP) input.solid input.dat	
	$(CP) corr_energy.txt Solid_result/corr_energy.txt
	$(CP) corr_pressure.txt Solid_result/corr_pressure.txt
	$(CP) output.gave.0 Solid_result/output.gave.0 
	$(CP) output.pres.0 Solid_result/output.press.0
	$(CP) output.epot.0 Solid_result/output.epot.0
Liquid: input.liquid input.dat# Copio i dati del liquido in quelli di input e copio gli output nella relativa cartella
	$(CP) input.liquid input.dat
	$(CP) corr_energy.txt Liquid_result/corr_energy.txt
	$(CP) corr_pressure.txt Liquid_result/corr_pressure.txt
	$(CP) output.gave.0 Liquid_result/output.gave.0 
	$(CP) output.pres.0 Liquid_result/output.press.0
	$(CP) output.epot.0 Liquid_result/output.epot.0
Gas: input.gas input.dat # Copio i dati del gas in quelli di input e copio gli output nella relativa cartella
	$(CP) input.liquid input.dat
	$(CP) input.gas input.dat
	$(CP) corr_energy.txt Gas_result/corr_energy.txt
	$(CP) corr_pressure.txt Gas_result/corr_pressure.txt
	$(CP) output.gave.0 Gas_result/output.gave.0 
	$(CP) output.pres.0 Gas_result/output.press.0
	$(CP) output.epot.0 Gas_result/output.epot.0

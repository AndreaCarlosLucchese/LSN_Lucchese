LIBS:=`root-config --libs`
INCS:=`root-config --cflags`
CP =cp
%.o : %.cpp
	g++ -Wall -c $< ${INCS}

qmc1d: qmc1d.o
	g++ -O3 -Wall -o $@ $^ ${LIBS}

clean:
	rm *.o qmc1d potential.dat kinetic.dat probability.dat

pigs: input.pigs input.dat
	$(CP) input.pigs input.dat

pimc: input.pimc input.dat
	$(CP) input.pimc input.dat

constant:
	$(CP) potential.dat constant_psi/potential.dat
	$(CP) kinetic.dat constant_psi/kinetic.dat
	$(CP) probability.dat constant_psi/probability.dat

GS:
	$(CP) potential.dat GS_psi/potential.dat
	$(CP) kinetic.dat GS_psi/kinetic.dat
	$(CP) probability.dat GS_psi/probability.dat

Pimc:
	$(CP) potential.dat Pimc_psi/potential.dat
	$(CP) kinetic.dat Pimc_psi/kinetic.dat
	$(CP) probability.dat Pimc_psi/probability.dat

TR:
	$(CP) potential.dat TR_psi/potential.dat
	$(CP) kinetic.dat TR_psi/kinetic.dat
	$(CP) probability.dat TR_psi/probability.dat

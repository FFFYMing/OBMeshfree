CC=gcc
CPP=g++
CFLAGS=-g -Wall

PDDIR=./bond-based-pd
NLDIFFDIR=./nonlocal_diffusion
KWDIR=./Kalthoff-Winkler-test

LAPACKFLAG=-llapack 

PD:
	$(CPP) $(PDDIR)/PMB_2Dweight.cpp $(LAPACKFLAG) -o PMB_2d.ex 
Nldiff:
	$(CPP) $(NLDIFFDIR)/nonlocaldiff_static.cpp $(LAPACKFLAG) -o nldiff.ex
Nldiffnl:
	$(CPP) $(NLDIFFDIR)/nonlocaldiff_static_nonlocal.cpp $(LAPACKFLAG) -o nldiffnl.ex
Nldiffd:
	$(CPP) $(NLDIFFDIR)/nonlocaldiff.cpp $(LAPACKFLAG) -o nldiffd.ex
KW:
	$(CPP) $(KWDIR)/KW_2Dweight_dynamic.cpp $(LAPACKFLAG) -o KW.ex
clean:
	rm -f *.ex


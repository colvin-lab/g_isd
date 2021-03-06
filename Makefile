### This works for a default install of gromacs from source.
### Gromacs library linking...
GMXINC=/usr/local/gromacs-4.6.5/include
GMXLIB=/usr/local/gromacs-4.6.5/lib
OPTIONS=-g -O2 -fopenmp -I$(GMXINC) -L$(GMXLIB)
LIBLINKS=-lgmx -lgmxana -lm
INSTALL_BINS=g_isd g_isdcalc g_isdcmds g_isdorder
INSTALL_PATH=$(HOME)/bin

###
### Basic compiler setup
###
CC=gcc
CPP=g++
F77=gfortran

# Targets
all : g_isd g_isdcalc g_isdcmds g_isdorder

# Executables
g_isd : g_isd.o gmx_isd.o libisdm.o libesa.o
	$(CC) $(OPTIONS) $^ -o $@ $(LIBLINKS)

g_isdcalc : g_isdcalc.o gmx_isdcalc.o libisdm.o libesa.o
	$(CC) $(OPTIONS) $^ -o $@ $(LIBLINKS)

g_isdcmds : g_isdcmds.o gmx_isdcmds.o libisdm.o libesa.o
	$(CC) $(OPTIONS) $^ -o $@ $(LIBLINKS)

g_isdorder : g_isdorder.o gmx_isdorder.o libisdm.o libesa.o
	$(CC) $(OPTIONS) $^ -o $@ $(LIBLINKS)

# Libraries

# Objects
g_isd.o : g_isd.c
	$(CC) $(OPTIONS) -c g_isd.c -o $@

gmx_isd.o : gmx_isd.c
	$(CC) $(OPTIONS) -c $^ -o $@

g_isdcalc.o : g_isdcalc.c
	$(CC) $(OPTIONS) -c g_isdcalc.c -o $@

gmx_isdcalc.o : gmx_isdcalc.c
	$(CC) $(OPTIONS) -c $^ -o $@

g_isdcmds.o : g_isdcmds.c
	$(CC) $(OPTIONS) -c g_isdcmds.c -o $@

gmx_isdcmds.o : gmx_isdcmds.c
	$(CC) $(OPTIONS) -c $^ -o $@

g_isdorder.o : g_isdorder.c
	$(CC) $(OPTIONS) -c g_isdorder.c -o $@

gmx_isdorder.o : gmx_isdorder.c
	$(CC) $(OPTIONS) -c $^ -o $@

libisdm.o : libisdm.c
	$(CC) $(OPTIONS) -c $^ -o $@

libesa.o : libesa.c
	$(CC) $(OPTIONS) -c $^ -o $@

clean :
	rm *.o

dist-clean : clean
	rm g_isd g_isdcalc g_isdcmds g_isdorder

install :
	mkdir -p $(INSTALL_PATH)
	cp $(INSTALL_BINS) $(INSTALL_PATH)

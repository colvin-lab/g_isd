### This works for a default install of gromacs from source.
### Gromacs library linking...
GMXINC=/usr/local/gromacs/include
GMXLIB=/usr/local/gromacs/lib
OPTIONS=-g -I$(GMXINC) -L$(GMXLIB)

###
### Basic compiler setup
###
CC=gcc
CPP=g++
F77=gfortran

# Targets
all : g_isdmap g_isdcalc g_isddecorr

# Executables
g_isdmap : g_isdmap.o gmx_isdmap.o libisdm.o libesa.o
	$(CC) $(OPTIONS) $^ -o $@ -lgmx -lgmxana

g_isdcalc : g_isdcalc.o gmx_isdcalc.o libisdm.o libesa.o
	$(CC) $(OPTIONS) $^ -o $@ -lgmx -lgmxana

g_isddecorr : g_isddecorr.o gmx_isddecorr.o libisdm.o libesa.o
	$(CC) $(OPTIONS) $^ -o $@ -lgmx -lgmxana

# Libraries

# Objects
g_isdmap.o : g_isdmap.c
	$(CC) $(OPTIONS) -c g_isdmap.c -o $@

gmx_isdmap.o : gmx_isdmap.c
	$(CC) $(OPTIONS) -c $^ -o $@

g_isdcalc.o : g_isdcalc.c
	$(CC) $(OPTIONS) -c g_isdcalc.c -o $@

gmx_isdcalc.o : gmx_isdcalc.c
	$(CC) $(OPTIONS) -c $^ -o $@

g_isddecorr.o : g_isddecorr.c
	$(CC) $(OPTIONS) -c g_isddecorr.c -o $@

gmx_isddecorr.o : gmx_isddecorr.c
	$(CC) $(OPTIONS) -c $^ -o $@

libisdm.o : libisdm.c
	$(CC) $(OPTIONS) -c $^ -o $@

libesa.o : libesa.c
	$(CC) $(OPTIONS) -c $^ -o $@

clean :
	rm *.o

dist-clean : clean
	rm g_isdmap g_isdcalc g_isddecorr

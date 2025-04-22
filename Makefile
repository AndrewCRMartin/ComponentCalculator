CC     = cc
COPT   = -ansi -pedantic -Wall -I$(HOME)/include
LOPT   = -L$(HOME)/lib
LIBS   = -lm
EXE    = compcalc
OFILES = compcalc.o

$(EXE) : $(OFILES)
	$(CC) $(LOPT) -o $(EXE) $(OFILES) $(LIBS)

.c.o :
	$(CC) $(COPT) -c -o $@ $<



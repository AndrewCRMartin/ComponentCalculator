CC     = cc
COPT   = -ansi -pedantic -Wall -I$(HOME)/include -g
LOPT   = -L$(HOME)/lib -g
LIBS   = -lm
EXE    = compcalc
OFILES = main.o compcalc.o
HFILES = compcalc.h eseries.h

$(EXE) : $(OFILES)
	$(CC) $(LOPT) -o $(EXE) $(OFILES) $(LIBS)

main.o : main.c $(HFILES)
	$(CC) $(COPT) -c -o $@ $<

test.o : test.c $(HFILES)
	$(CC) $(COPT) -c -o $@ $<

compcalc.o : compcalc.c $(HFILES)
	$(CC) $(COPT) -c -o $@ $<

.c.o :
	$(CC) $(COPT) -c -o $@ $<

test : test.o compcalc.o
	$(CC) $(LOPT) -o test test.o compcalc.o $(LIBS)

clean :
	\rm $(OFILES) test.o

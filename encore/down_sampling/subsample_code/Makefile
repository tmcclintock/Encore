#makefile
CC          =  gcc
OPTIMIZE    =  -g -O3 #-Wall -wd981 #-wd1419 -wd810
OPENMPFLAG  =  
GSLI        =  -I/home/tmcclintock/code/gsl/include
GSLL        =  -L/home/tmcclintock/code/gsl/lib -lgsl -lgslcblas
EXTRACFLAGS =  -Wall -W -Wmissing-prototypes -Wstrict-prototypes -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align \
	-Wwrite-strings -Wnested-externs -fshort-enums -fno-common -Dinline=
EXTRACLIB   =

CLINK=$(CC)
CFLAGS=$(OPTIMIZE) $(GSLI) $(EXTRACFLAGS) $(OPTS) 
CLIB=$(EXTRACLIB) $(GSLL) -lm

OBJS = io.o io_LGADGET.o 

EXEC = subsamp_parts
TEST =
all: $(EXEC) 
test: $(TEST) clean

OBJS1=$(OBJS) main.o
$(EXEC): $(OBJS1)
	$(CLINK) $(CFLAGS) -o $@ $(OBJS1) $(CLIB)

$(OBJS1): Makefile allheader.h 

.PHONY : clean
clean: 
	rm -f *.o

.PHONY : spotless
spotless: 
	rm -f *.o $(EXEC) $(TEST)

.PHONY : pristine
pristine: 
	rm -f *.o $(EXEC) $(TEST) *~


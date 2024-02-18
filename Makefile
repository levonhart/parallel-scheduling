SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic
CPLEXROOT := /opt/ibm/ILOG/CPLEX_Studio2211
CPLEXDIR := $(CPLEXROOT)/cplex
CPLEXINCDIR := $(CPLEXDIR)/include
CONCERTDIR := $(CPLEXROOT)/concert
CONCERTINCDIR := $(CONCERTDIR)/include
CPLEXLIB      = cplex$(dynamic:yes=2211)

CC := clang
CCC := clang++
DFLAGS := -g -gdwarf
CCFLAGS := $(DFLAGS) -I../include
LNDIRS := -L$(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT) -L$(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
LDFLAGS := -L../libs -lm -lpthread -ldl -lconcert -lilocplex -l$(CPLEXLIB)
CPLEXFLAGS := -I$(CPLEXINCDIR) -I$(CONCERTINCDIR)

TARGETS:= cutting
MAIN   := cutting.o
OBJ    := cutting.o
DEPS   :=

.PHONY: all, clean

all: $(TARGETS)

$(OBJ): %.o : %.cpp $(DEPS)
		$(CCC) -c -o $@ $< $(CCFLAGS) $(CPLEXFLAGS)

$(TARGETS): % : $(filter-out $(MAIN), $(OBJ)) %.o
		$(CCC) -o $@ $(LIBS) $^ $(CCFLAGS) $(LNDIRS) $(LDFLAGS)

clean:
		rm -f $(TARGETS) $(OBJ)


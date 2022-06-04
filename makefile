SDIR = src

IDIR = include
CCCMD = gcc
CFLAGS = -I$(IDIR) -Wall

debug: CC = $(CCCMD) -D__GRAPH_C_DEBUG_ -D__GRAPH_C_DETECTION_DEBUG_
debug: BDIR = debug

release: CC = $(CCCMD) -O2
release: BDIR = build

ODIR=.obj
LDIR=lib

LIBS = -lm

_DEPS = graph.h fluid_mechanics.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = main.o graph.o fluid_mechanics.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: $(SDIR)/%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

release: $(OBJ)
	mkdir -p $(ODIR)
	mkdir -p $(BDIR)
	$(CC) -o $(BDIR)/LeakDetection $^ $(CFLAGS) $(LIBS)

debug: $(OBJ)
	mkdir -p $(ODIR)
	mkdir -p $(BDIR)
	$(CC) -o $(BDIR)/LeakDetection $^ $(CFLAGS) $(LIBS)

.PHONY: clean
clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~

.PHONY: all
all: main clean

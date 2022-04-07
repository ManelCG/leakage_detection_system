BDIR = build
SDIR = src

IDIR = include
CC = gcc
CFLAGS = -I$(IDIR) -Wall

ODIR=.obj
LDIR=lib

LIBS = -lm

_DEPS = graph.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = main.o graph.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: $(SDIR)/%.c $(DEPS)
	mkdir -p $(ODIR)
	$(CC) -c -o $@ $< $(CFLAGS)

main: $(OBJ)
	mkdir -p $(BDIR)
	$(CC) -o $(BDIR)/$@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean
clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~

.PHONY: all
all: main clean
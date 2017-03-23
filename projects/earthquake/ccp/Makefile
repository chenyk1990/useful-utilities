GMT_INC = -I/opt/local/lib/gmt4/include -I/opt/local/include
GMT_LIBS = -L/opt/local/lib/gmt4/lib -lgmt -lpsl -L/opt/local/lib -lnetcdf -lm -s

CFLAGS=-O ${GMT_INC}

SUBS = fft.o Complex.o Intp2D.o sacio.o

all: ccpStack3D 3DSlice

ccpStack3D: ccpStack3D.o $(SUBS)
	$(LINK.c) -o $@ $@.o $(SUBS) -lm

3DSlice: 3DSlice.o $(SUBS)
	$(LINK.c) -o $@ $@.o $(SUBS) ${GMT_LIBS}

clean:
	rm -f *.o ccpStack3D 3DSlice

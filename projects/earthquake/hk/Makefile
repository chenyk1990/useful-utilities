CFLAGS=-O ${GMT_INC}

SUBS = sacio.o fft.o Complex.o

all: rcvFn decon iter_decon k_stack grdmin

rcvFn: rcvFn.o knt.o matrix.o $(SUBS)
	$(LINK.c) -o $@ $@.o knt.o matrix.o $(SUBS) -lm

decon: decon.o decon_sub.o $(SUBS)
	$(LINK.c) -o $@ $@.o decon_sub.o $(SUBS) -lm

iter_decon: iter_decon.o decon_sub.o $(SUBS)
	$(LINK.c) -o $@ $@.o decon_sub.o $(SUBS) -lm


k_stack: k_stack.o $(SUBS)
	$(LINK.c) -o $@ $@.o $(SUBS) ${GMT_LIBS}

grdmin: grdmin.o grid2d.o
	$(LINK.c) -o $@ $@.o grid2d.o ${GMT_LIBS}

clean:
	rm -f *.o rcvFn decon iter_decon k_stack grdmin

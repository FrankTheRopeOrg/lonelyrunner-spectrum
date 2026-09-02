CC=gcc
CFLAGS=-O2 -fopenmp
all: mltest fastrunner fastrunner5 prejump ucscan ucdump ucpoints survey
%: %.c lrk.h
	$(CC) $(CFLAGS) -o $@ $<
clean:
	rm -f mltest fastrunner fastrunner5 prejump ucscan ucdump ucpoints survey

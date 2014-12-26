
This is a program which creates a map of the effects of background selection (BS) given a genetic map, locations of conserved elements and selection parameters.
The code was originally written by Graham McVicker (see McVicker et al. 2009), wrapped by Eyal Elyashiv and adapted to unix by Guy Amster.


Dependencies:

Windows:
1. gsl (libgsl.a)
2. glib 2.26.1 (probably supports 2.0 also) (glib-2.0.lib)
3. zlib 1.2.5 (zlib.lib)

Unix:
1. gsl 1.6 (with the supporting CBLAS library, -lgslcblas)
2. glib 2.0
3. zlib
4. math library (-lm).

NOTE: the Makefile is NOT generic, see these 2 lines - 
MY_CFLAGS = -I/ifs/data/c2b2/gs_lab/shared/software/gsl-1.16/include/gsl/ $(shell pkg-config --cflags glib-2.0) 
MY_LIBS   = -L/ifs/data/c2b2/gs_lab/shared/software/gsl-1.16/lib/ -lgsl -lgslcblas -lm -lz $(shell pkg-config --libs glib-2.0) 


Eyal Elyashiv
August 2014
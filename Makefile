LIBSRC := $(wildcard src/*.f)
LIBOBJ := $(patsubst %.f,%.o, $(LIBSRC))
BINSRC := $(wildcard bin/*.f)
F77    := gfortran
LIBSTA := libbifpack.a
LIBDYN := libbifpack.so
LIBDIR := lib
MKDIR  := mkdir -p

.SUFFIXES: .o .f

.f.o:
	$(F77) -c -o $@ $<

static: $(LIBOBJ) $(LIBDIR)
	ar cr $(LIBDIR)/$(LIBSTA) $(LIBOBJ)

$(LIBDIR):
	$(MKDIR) $(LIBDIR)

clean:
	$(RM) $(LIBOBJ)

distclean: clean
	$(RM) -r $(LIBDIR)


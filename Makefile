LIBSRC := $(wildcard src/*.f)
BIFOBJ := $(patsubst %.f,%.o, $(LIBSRC))
BINSRC := $(wildcard bin/*.f)
BIFBIN := $(addprefix bin/bifpack-, $(filter-out dataman, $(notdir $(patsubst %.f,%, $(BINSRC)))))
SOLSRC := $(filter-out solvers/shoot.f, $(wildcard solvers/*.f))
SOLOBJ := $(patsubst %.f,%.o, $(SOLSRC))
EXASRC := $(wildcard examples/exa*.f)
EXABIN := $(patsubst %.f,%, $(EXASRC))
EXBSRC := $(wildcard examples/exb*.f)
EXBBIN := $(patsubst %.f,%, $(EXBSRC))
EXDSRC := $(wildcard examples/exd*.f)
EXDBIN := $(patsubst %.f,%, $(EXDSRC))
F77    := gfortran -std=legacy
LIBBIF := bifpack
LIBSOL := bifsolv
LIBDIR := lib
MKDIR  := mkdir -p

.SUFFIXES: .o .f

.f.o:
	$(F77) -c -o $@ $<

bin/bifpack-%: bin/%.f
	$(F77) -o $@ -L ./$(LIBDIR) $<

examples/exa%: examples/exa%.f static
	$(F77) -o $@ $< examples/maina.f -L./$(LIBDIR) -l$(LIBBIF) -l$(LIBSOL)

examples/exb%: examples/exb%.f static
	$(F77) -o $@ $< examples/mainb.f -L./$(LIBDIR) -l$(LIBBIF) -l$(LIBSOL)

examples/exd%: examples/exd%.f static
	$(F77) -o $@ $< examples/maind.f -L./$(LIBDIR) -l$(LIBBIF) -l$(LIBSOL)

static: $(BIFOBJ) $(SOLOBJ) $(LIBDIR)
	ar cr $(LIBDIR)/lib$(LIBBIF).a $(BIFOBJ)
	ar cr $(LIBDIR)/lib$(LIBSOL).a $(SOLOBJ)

.PHONY: bin
bin: static $(BIFBIN)

$(LIBDIR):
	$(MKDIR) $(LIBDIR)

clean:
	$(RM) $(BIFOBJ) $(SOLOBJ)

distclean: clean
	$(RM) $(BIFBIN) 
	$(RM) $(EXABIN) 
	$(RM) $(EXBBIN) 
	$(RM) $(EXDBIN)
	$(RM) -r $(LIBDIR) 


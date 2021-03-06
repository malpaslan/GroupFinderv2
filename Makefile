# Make file for GroupFinderv2.

# Basic definitions.
CC = gcc
ODIR = obj
IDIR = include
LDIR = /home/users/ma5046/libC_main
HOME = /home/users/ma5046/GroupFinderv2
CFLAGS = -fopenmp -L$(LDIR) -I$(IDIR) -nostartfiles

# Libraries, objects, and dependencies.
LIBS = -lC_main -lm
_OBJ = spline.o splint.o qromo.o midpnt.o polint.o sort2.o sort3.o sham.o zbrent.o trapzd.o qtrap.o scatter.o gammln.o fiber_corrected_galaxy_property.o
OBJ = $(patsubst %, $(ODIR)/%,$(_OBJ))

_DEPS = kdtree.h header.h nrutil.h
DEPS = $(patsubst %, $(IDIR)/%,$(_DEPS))

# Build object files.
$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

# Compile groupfinding codes.

GroupFinderv2Main: $(OBJ)
	$(CC) -o $@ $(CFLAGS) $(LIBS)
	cp -f $@ $(HOME)/$@

GroupFinderv2Mocks: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)
	cp -f $@ $(HOME)/$@

GroupFinderv2Mocks_ColWeight: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)
	cp -f $@ $(HOME)/$@

GroupFinderv2Mocks_ColWeight_OMP: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)
	cp -f $@ $(HOME)/$@

GroupFinderv2Box_ColWeight_OMP: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)
	cp -f $@ $(HOME)/$@

GroupFinderv2_2MRS: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)
	cp -f $@ $(HOME)/$@

# Cleanup.
.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 


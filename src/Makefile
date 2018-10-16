# Make file for GroupFinderv2.

# Basic definitions.
CC = gcc
LDIR = /home/users/ma5046/libC_main
HOME = /home/users/ma5046/GroupFinderv2
ODIR = $(HOME)/src/obj
IDIR = $(HOME)/src/include

CFLAGS = -fopenmp -L$(LDIR) -I$(IDIR)

# Libraries, objects, and dependencies.
LIBS = -lC_main -lm
_OBJ = spline.o splint.o qromo.o midpnt.o polint.o sort2.o sort3.o sham.o zbrent.o trapzd.o qtrap.o scatter.o gammln.o fiber_corrected_galaxy_property.o kdtree.o header.o 
OBJ = $(patsubst %, $(ODIR)/%,$(_OBJ))

_DEPS = kdtree.h header.h nrutil.h
DEPS = $(patsubst %, $(IDIR)/%,$(_DEPS))

# Executable list.

EXECS = GroupFinderv2Main GroupFinderv2Mocks GroupFinderv2Mocks_ColWeight GroupFinderv2Mocks_ColWeight_OMP GroupFinderv2Box_ColWeight_OMP GroupFinderv2_2MRS

# Build object files.
$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

# Compile groupfinding codes.
.PHONY: all

all: $(EXECS) 

GroupFinderv2Main: GroupFinderv2Main.c $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)
	cp -f $@ $(HOME)/$@

GroupFinderv2Mocks: GroupFinderv2Mocks.c $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)
	cp -f $@ $(HOME)/$@

GroupFinderv2Mocks_ColWeight: GroupFinderv2Mocks_ColWeight.c $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)
	cp -f $@ $(HOME)/$@

GroupFinderv2Mocks_ColWeight_OMP: GroupFinderv2Mocks_ColWeight_OMP.c $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)
	cp -f $@ $(HOME)/$@

GroupFinderv2Box_ColWeight_OMP: GroupFinderv2Box_ColWeight_OMP.c $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)
	cp -f $@ $(HOME)/$@

GroupFinderv2_2MRS: GroupFinderv2_2MRS.c $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)
	cp -f $@ $(HOME)/$@

# Cleanup.
.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 


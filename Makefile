HOME = /home/users/ma5046/
LIB = -L${HOME}libC_main -lC_main -lm
GFDIR = ${HOME}/GroupFinderv2/

CC = gcc
CFLAGS = -O2

OBJSGF = fiber_corrected_galaxy_property.o sham.o gammln.o sort2.o kdtree.o sort3.o midpnt.o spline.o polint.o splint.o qromo.o trapzd.o qtrap.o zbrent.o scatter.o
GroupFinderv2Mocks_ColWeight: $(OBJSGF)
	$(CC) -o gfv2mockCol $(OBJSGF) $(LIB)

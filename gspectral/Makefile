########################################################################
# GOOPS - Generic Object Oriented Pseudo Spectra
# gspectral subdir	
# Manuel Baptista (2004) (2005)
########################################################################

CXX = g++ -O2 -mtune=pentium-m

########################################################################

INCLUDE = -I/home/mbaptist/work/codes/devel/cat -I..

OBJECTS = spectral_fourier.o spectral_fourier_layer.o

########################################################################

all: $(OBJECTS)

%.o: %.C *.h
	$(CXX) $(INCLUDE) -c $<

clean:
	@rm -rfv *~ *.o


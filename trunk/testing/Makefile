########################################################################
# GOOPS - Generic Object Oriented Pseudo Spectra
# testing
# Manuel Baptista (2004) (2005)
########################################################################

CXX = g++ -O3 -funroll-loops -fexpensive-optimizations -g

INCLUDE = -I/home/mbaptist/work/codes/devel/cat -I..

########################################################################


all: testgffti

testgffti: testgffti.o
	$(CXX) -o testgffti testgffti.o -L.. -lgoops -lfftw3

%.o: %.C
	$(CXX) $(INCLUDE) -c $<

clean:
	@rm -rfv *~ *.o testgffti



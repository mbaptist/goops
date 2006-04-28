########################################################################
# GOOPS - Generic Object Oriented Pseudo Spectra
# gffti subdir	
# Manuel Baptista (2004) (2005)
########################################################################

CXX = g++ -O3 

########################################################################

INCLUDE = -I/home/mbaptist/work/codes/devel/cat 

OBJECTS = s_rfft_h_z.o v_rfft_h_z.o

########################################################################

all: $(OBJECTS)

%.o: %.C *.h
	$(CXX) $(INCLUDE) -c $<

clean:
	@rm -rfv *~ *.o




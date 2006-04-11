########################################################################
# GOOPS - Generic Object Oriented Pseudo Spectra
# gffti subdir	
# Manuel Baptista (2004) (2005)
########################################################################

CXX = g++ -g -O2

########################################################################

INCLUDE = -I/home/mbaptist/work/codes/devel/cat 

OBJECTS = v_rfft_h_z_ssc.o v_rfft_h_z_ccs.o

########################################################################

all: $(OBJECTS)

%.o: %.C *.h
	$(CXX) $(INCLUDE) -c $<

clean:
	@rm -rfv *~ *.o




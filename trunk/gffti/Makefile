########################################################################
# GOOPS - Generic Object Oriented Pseudo Spectra
# gffti subdir	
# Manuel Baptista (2004) (2005)
########################################################################

CXX = g++ -O3 -m32 -msse2 -mfpmath=sse

########################################################################

INCLUDE = -I/home/mbaptist/work/codes/devel/cat 

OBJECTS = s_cosfft_1d.o s_sinfft_1d.o v_rfft_h_z_ssc.o v_rfft_h_z_ccs.o

########################################################################

all:
	g++ -I../../cat/ testgffti.C -o testgffti

#$(OBJECTS)

#%.o: %.C *.h
#	$(CXX) $(INCLUDE) -c $<

clean:
	@rm -rfv *~ *.o




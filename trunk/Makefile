########################################################################
# GOOPS - Generic Object Oriented Pseudo Spectra
# Manuel Baptista (2004) (2005)
########################################################################

CXX = g++ -O3 -funroll-loops -fexpensive-optimizations -g

INCLUDE = -I/home/mbaptist/work/codes/devel/cat 

LIB = -lpthread 

########################################################################

GFFTI_OBJECTS = \
  ./gffti/s_rfft_h_z.o \
	./gffti/v_rfft_h_z.o \

GSPECTRAL_OBJECTS= \
        ./gspectral/spectral_fourier.o \
        ./gspectral/spectral_fourier_layer.o

OBJECTS = $(GFFTI_OBJECTS) $(GSPECTRAL_OBJECTS)

GFFTI_INSTALL_HEADERS = \
	./gffti/gffti.h \
	./gffti/s_fft.h \
	./gffti/s_rfft.h \
	./gffti/s_sinfft_1d.h \
	./gffti/s_cosfft_1d.h \
	./gffti/s_fft_h_z.h \
	./gffti/s_rfft_h_z.h \
	./gffti/v_rfft.h \
	./gffti/v_rfft_h_z_ccs.h \
	./gffti/v_rfft_h_z_ssc.h

GSPECTRAL_INSTALL_HEADERS = \
	./gspectral/gspectral.h \
	./gspectral/spectral_fourier.h \
	./gspectral/spectral_fourier_layer.h
		
INSTALL_HEADERS = $(GFFTI_INSTALL_HEADERS) $(GSPECTRAL_INSTALL_HEADERS)

INSTALL_ROOT = /usr/local

########################################################################

all: libgoops 
#test

libgoops:
	$(MAKE) -C ./gffti
	$(MAKE) -C ./gspectral
	$(CXX) -shared -o libgoops.so $(OBJECTS) $(LIB)
	ar rcs libgoops.a $(OBJECTS)

test:
	$(MAKE) -C ./testing 

clean:
	@rm -rfv *~ *.o *.a *.so
	$(MAKE) -C ./gffti clean 
	$(MAKE) -C ./gspectral clean 
	$(MAKE) -C ./testing clean 

install:
	install -c -m 644 $(GFFTI_INSTALL_HEADERS) \
		$(INSTALL_ROOT)/include/goops/gffti
	install -c -m 644 $(GSPECTRAL_INSTALL_HEADERS) \
		$(INSTALL_ROOT)/include/goops/gspectral
	install -c -m 644 goops.h $(INSTALL_ROOT)/include/goops
	install -c -m 755 libgoops.so $(INSTALL_ROOT)/lib 
	install -c -m 644 libgoops.a $(INSTALL_ROOT)/lib
	@/sbin/ldconfig

uninstall:
	@rm -vrf $(INSTALL_ROOT)/lib/libgoops* \
		$(INSTALL_ROOT)/include/goops
	@/sbin/ldconfig


########################################################################
# GOOPS - Generic Object Oriented Pseudo Spectra
# Manuel Baptista (2004) (2005)
########################################################################

CXX = g++ -O3 -funroll-loops -fexpensive-optimizations -g 

INCLUDE = -I/home/mbaptist/work/codes/devel/cat 

LIB = -lpthread 

########################################################################

FFT_OBJECTS = \
  ./fft/s_rfft_h_z.o \
	./fft/v_rfft_h_z.o

SPECTRAL_OBJECTS= \
        ./spectral/spectral_fourier_base.o \
        ./spectral/spectral_fourier.o \
        ./spectral/spectral_fourier_layer.o

OBJECTS = $(FFT_OBJECTS) $(SPECTRAL_OBJECTS)

FFT_INSTALL_HEADERS = \
	./fft/fft.h \
	./fft/s_fft.h \
	./fft/s_rfft.h \
	./fft/s_sinfft_1d.h \
	./fft/s_cosfft_1d.h \
	./fft/s_fft_h_z.h \
	./fft/s_rfft_h_z.h \
	./fft/v_rfft.h \
	./fft/v_rfft_h_z_ccs.h \
	./fft/v_rfft_h_z_ssc.h

SPECTRAL_INSTALL_HEADERS = \
	./spectral/spectral.h \
  ./spectral/spectral_fourier_base.h \
	./spectral/spectral_fourier.h \
	./spectral/spectral_fourier_layer.h
		
INSTALL_HEADERS = $(FFT_INSTALL_HEADERS) $(SPECTRAL_INSTALL_HEADERS)

INSTALL_ROOT = /usr/local

########################################################################

all: libgoops 
#test

libgoops:
	$(MAKE) -C ./fft
	$(MAKE) -C ./spectral
	$(CXX) -shared -o libgoops.so $(OBJECTS) $(LIB)
	ar rcs libgoops.a $(OBJECTS)

test:
	$(MAKE) -C ./testing 

clean:
	@rm -rfv *~ *.o *.a *.so
	$(MAKE) -C ./fft clean
	$(MAKE) -C ./spectral clean
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


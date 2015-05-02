project(libgoops)
cmake_minimum_required(VERSION 2.8)


set(CAT_INSTALL_DIR "${CMAKE_SOURCE_DIR}/../cat/install")
include_directories("${CAT_INSTALL_DIR}/include/cat")
link_directories("${CAT_INSTALL_DIR}/lib")


set(GOOPS_DIR "${CMAKE_SOURCE_DIR}/goops")

include_directories("${GOOPS_DIR}")
include_directories("${GOOPS_DIR}/fft")
include_directories("${GOOPS_DIR}/spectral")


set(GOOPS_HEADERS

"${GOOPS_DIR}/goops.h"
"${GOOPS_DIR}/goops_types.h"

)


set(FFT_HEADERS

"${GOOPS_DIR}/fft/fftclass.h"
"${GOOPS_DIR}/fft/fft.h"
"${GOOPS_DIR}/fft/fft_types.h"
"${GOOPS_DIR}/fft/plan.h"
"${GOOPS_DIR}/fft/s_fft_h_z.h"
"${GOOPS_DIR}/fft/s_rfft_h_z.h"
"${GOOPS_DIR}/fft/v_rfft.h"
"${GOOPS_DIR}/fft/v_rfft_h_z.h"

)

set(SPECTRAL_HEADERS

"${GOOPS_DIR}/spectral/globals.h"  
"${GOOPS_DIR}/spectral/spectral_fourier_base.h"  
"${GOOPS_DIR}/spectral/spectral_fourier.h"  
"${GOOPS_DIR}/spectral/spectral_fourier_layer.h"  
"${GOOPS_DIR}/spectral/spectral.h"

)

set(FFT_TEMPLATE_IMPLEMENTATIONS

"${GOOPS_DIR}/fft/fftclass.C"
"${GOOPS_DIR}/fft/plan.C"

)

set(FFT_SOURCES

"${GOOPS_DIR}/fft/s_rfft_h_z.C"
"${GOOPS_DIR}/fft/v_rfft_h_z.C"

)

set(SPECTRAL_SOURCES

"${GOOPS_DIR}/spectral/spectral_fourier_base.C"  
"${GOOPS_DIR}/spectral/spectral_fourier.C"  
"${GOOPS_DIR}/spectral/spectral_fourier_layer.C" 

)


set(LIB_TARGET goops)
add_library(${LIB_TARGET} SHARED ${FFT_SOURCES} ${SPECTRAL_SOURCES})
target_link_libraries(goops cat)

install (TARGETS ${LIB_TARGET} DESTINATION lib)
install (FILES ${GOOPS_HEADERS} DESTINATION include/${LIB_TARGET})
install (FILES ${FFT_HEADERS} DESTINATION include/${LIB_TARGET}/fft)
install (FILES ${SPECTRAL_HEADERS} DESTINATION include/${LIB_TARGET}/spectral)
install (FILES ${FFT_TEMPLATE_IMPLEMENTATIONS} DESTINATION include/${LIB_TARGET}/fft)


set(TESTGOOPS_TARGET testgoops)

enable_testing()

find_package (Cat)
if (CAT_FOUND)
  include_directories(${CAT_INCLUDE_DIR})
  set(LIBS ${LIBS} ${CAT_LIBRARIES})
endif (CAT_FOUND)

find_package (FFTW)
if (FFTW_FOUND)
  include_directories(${FFTW_INCLUDE_DIR})
  set(LIBS ${LIBS} ${FFTW_LIBRARIES})
endif (FFTW_FOUND)

set(TESTGOOPS_DIR "${CMAKE_SOURCE_DIR}/test")
include_directories("${TESTGOOPS_DIR}")

set(TESTGOOPS_SOURCES 

#"${TESTGOOPS_DIR}/testgffti.C"
"${TESTGOOPS_DIR}/testsffthz.C"

)
add_executable(${TESTGOOPS_TARGET} ${TESTGOOPS_SOURCES})
target_link_libraries(${TESTGOOPS_TARGET} ${LIBS} goops)
add_dependencies(${TESTGOOPS_TARGET} goops)

add_test(${TESTGOOPS_TARGET} ${CMAKE_CURRENT_BINARY_DIR}/${TESTGOOPS_TARGET})




enable_testing()

set(LIBS goops cat)

set(CAT_INSTALL_DIR "${CMAKE_SOURCE_DIR}/../cat/install")
include_directories("${CAT_INSTALL_DIR}/include/cat")
link_directories("${CAT_INSTALL_DIR}/lib")

set(GOOPS_INSTALL_DIR "${CMAKE_SOURCE_DIR}/install")
include_directories("${GOOPS_INSTALL_DIR}/include/goops")
link_directories("${GOOPS_INSTALL_DIR}/lib")

set(GOOPS_DIR "${CMAKE_SOURCE_DIR}/goops")
include_directories("${GOOPS_DIR}")

set(TESTGOOPS_DIR "${CMAKE_SOURCE_DIR}/test")
include_directories("${TESTGOOPS_DIR}")

set(TESTGOOPS_SOURCES 

#"${TESTGOOPS_DIR}/testgffti.C"
"${TESTGOOPS_DIR}/testsffthz.C"

)

set(TESTGOOPS_TARGET testgoops)
add_executable(${TESTGOOPS_TARGET} ${TESTGOOPS_SOURCES})
target_link_libraries(${TESTGOOPS_TARGET} ${LIBS})
add_dependencies(${TESTGOOPS_TARGET} goops)

find_package (FFTW)
if (FFTW_FOUND)
  include_directories(${FFTW_INCLUDE_DIRS})
  target_link_libraries (testgoops ${FFTW_LIBRARIES})
endif (FFTW_FOUND)

add_test(${TESTGOOPS_TARGET} ${CMAKE_CURRENT_BINARY_DIR}/${TESTGOOPS_TARGET})


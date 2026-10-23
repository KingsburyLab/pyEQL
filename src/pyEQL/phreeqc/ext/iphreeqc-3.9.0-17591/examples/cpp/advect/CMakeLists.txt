# project
project(example_advect_cpp CXX)

configure_file(phreeqc.dat phreeqc.dat COPYONLY)
configure_file(ic ic COPYONLY)
configure_file(CMakeLists.txt.in CMakeLists.txt COPYONLY)

# files
SET(CPP_Advect_Files
${IPhreeqc_BINARY_DIR}/examples/cpp/advect/CMakeLists.txt
advect.cpp
ic
phreeqc.dat
README.txt
)

# src
SET(CPP_Advect_SRC
advect.cpp
)

# executable
add_executable(example_advect_cpp ${CPP_Advect_SRC})

# library dependencies
SET(EXTRA_LIBS ${EXTRA_LIBS} IPhreeqc)

# link
target_link_libraries(example_advect_cpp ${EXTRA_LIBS})

# install directory
SET(CPP_Advect_Dir ${EXAMPLES_DIR}/cpp/advect)

# install
install(FILES ${CPP_Advect_Files} DESTINATION ${CPP_Advect_Dir})

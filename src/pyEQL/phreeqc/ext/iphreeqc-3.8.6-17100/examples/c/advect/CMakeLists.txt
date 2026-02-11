# project
project(example_advect_c C)

configure_file(phreeqc.dat phreeqc.dat COPYONLY)
configure_file(ic ic COPYONLY)
configure_file(CMakeLists.txt.in CMakeLists.txt COPYONLY)

# files
SET(C_Advect_Files
${IPhreeqc_BINARY_DIR}/examples/c/advect/CMakeLists.txt
advect.c
ic
phreeqc.dat
README.txt
)

# src
SET(C_Advect_SRC
advect.c
)

# executable
add_executable(example_advect_c ${C_Advect_SRC})

# library dependencies
SET(EXTRA_LIBS ${EXTRA_LIBS} IPhreeqc)

# link
target_link_libraries(example_advect_c ${EXTRA_LIBS})

# install directory
SET(C_Advect_Dir ${EXAMPLES_DIR}/c/advect)

# install
install(FILES ${C_Advect_Files} DESTINATION ${C_Advect_Dir})

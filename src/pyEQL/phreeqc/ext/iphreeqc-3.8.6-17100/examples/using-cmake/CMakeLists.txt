configure_file(CMakeLists.txt.in CMakeLists.txt COPYONLY)

set(example_FILES
  ${IPhreeqc_BINARY_DIR}/examples/using-cmake/CMakeLists.txt
  ex2
  main.cpp
  phreeqc.dat
  README.txt
  )

# install directory
set(example_DIR ${EXAMPLES_DIR}/using-cmake)

if (IS_ABSOLUTE "${CMAKE_INSTALL_PREFIX}")
  set(FULL_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")
else()
  # this assumes that relative install prefixes are always relative to the build directory
  get_filename_component(FULL_INSTALL_PREFIX "${IPhreeqc_BINARY_DIR}/${CMAKE_INSTALL_PREFIX}" ABSOLUTE)
endif()
configure_file(post-install.cmake.in ${IPhreeqc_BINARY_DIR}/post-install.cmake @ONLY)

# install
install(FILES ${example_FILES} DESTINATION ${example_DIR})

# configure
install(CODE "execute_process(COMMAND ${CMAKE_CTEST_COMMAND} -S post-install.cmake)")

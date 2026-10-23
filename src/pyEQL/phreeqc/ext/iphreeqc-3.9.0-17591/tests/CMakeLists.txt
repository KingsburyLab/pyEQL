##
## Configure files
##

# tests/ex2
configure_file (
  "ex2.in"
  "ex2"
  )

# tests/phreeqc.dat
configure_file (
  "phreeqc.dat.in"
  "phreeqc.dat"
  )

##
## Test C
##

project(test_c C)

# source
SET(test_c_SOURCES
  test_c.c
)

# test executable
add_executable(test_c ${test_c_SOURCES})

# library dependencies
SET(EXTRA_LIBS ${EXTRA_LIBS} IPhreeqc)

# link
target_link_libraries(test_c ${EXTRA_LIBS})
if (MSVC AND BUILD_SHARED_LIBS AND BUILD_CLR_LIBS)
  # Setup references for /CLR
  set_target_properties(test_c PROPERTIES VS_DOTNET_REFERENCES "${CMAKE_CURRENT_SOURCE_DIR}/../src/phreeqcpp/ZedGraph.dll")
endif()


# test compile and run
add_test(TestCompileAndRunC test_c)

if (MSVC AND BUILD_SHARED_LIBS)
  # copy dll
  add_custom_command(TARGET test_c POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:IPhreeqc> $<TARGET_FILE_DIR:test_c>
  )
endif()

##
## Test CXX
##

project(test_cxx CXX)

# source
SET(test_cxx_SOURCES
  test_cxx.cxx
)

# test executable
add_executable(test_cxx ${test_cxx_SOURCES})

# library dependencies
SET(EXTRA_LIBS ${EXTRA_LIBS} IPhreeqc)

# link
target_link_libraries(test_cxx ${EXTRA_LIBS})
if (MSVC AND BUILD_SHARED_LIBS AND BUILD_CLR_LIBS)
  # Setup references for /CLR
  set_target_properties(test_cxx PROPERTIES VS_DOTNET_REFERENCES "${CMAKE_CURRENT_SOURCE_DIR}/../src/phreeqcpp/ZedGraph.dll")
endif()

# test compile and run
add_test(TestCompileAndRunCXX test_cxx)

if (MSVC AND BUILD_SHARED_LIBS)
  # copy dll
  add_custom_command(TARGET test_cxx POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:IPhreeqc> $<TARGET_FILE_DIR:test_cxx>
  )
endif()


##
## Test Fortran
##

if (IPHREEQC_FORTRAN_TESTING)

  project(test_f90 Fortran)

  # source
  SET(test_f90_SOURCES
    main.f90
    test_f90.F90
  )

  if (IPHREEQC_ENABLE_MODULE)
    SET(test_f90_SOURCES
      ${test_f90_SOURCES}
      ../src/IPhreeqc_interface.F90
    )
  endif()

  # place [*.F[90]] files into "Source Files"
  source_group("Source Files" FILES "test_f90.F90")
  if (IPHREEQC_ENABLE_MODULE)
    source_group("Source Files" FILES "../src/IPhreeqc_interface.F90")
  endif()

  # test executable
  add_executable(test_f90 ${test_f90_SOURCES})

  # library dependencies
  SET(EXTRA_LIBS ${EXTRA_LIBS} IPhreeqc)

  # link
  target_link_libraries(test_f90 ${EXTRA_LIBS})

  # test compile and run
  add_test(TestCompileAndRunFortran test_f90)

  if (MSVC AND BUILD_SHARED_LIBS)
    # copy dll
    add_custom_command(TARGET test_f90 POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:IPhreeqc> $<TARGET_FILE_DIR:test_f90>
    )
  endif()

endif()

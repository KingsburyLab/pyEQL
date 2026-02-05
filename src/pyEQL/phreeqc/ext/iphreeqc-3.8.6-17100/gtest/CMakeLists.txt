project(unit)

include(GoogleTest)

add_executable(TestCVar TestCVar.cpp)
target_link_libraries(TestCVar IPhreeqc gtest gtest_main gmock)
gtest_discover_tests(TestCVar)

add_executable(TestVar TestVar.cpp)
target_link_libraries(TestVar IPhreeqc gtest gtest_main gmock)
gtest_discover_tests(TestVar)

add_executable(TestIPhreeqcLib TestIPhreeqcLib.cpp FileTest.cpp)
target_link_libraries(TestIPhreeqcLib IPhreeqc gtest gtest_main gmock)
gtest_discover_tests(TestIPhreeqcLib)

add_executable(TestIPhreeqc TestIPhreeqc.cpp FileTest.cpp)
target_link_libraries(TestIPhreeqc IPhreeqc gtest gtest_main gmock)
if (MSVC AND BUILD_SHARED_LIBS AND BUILD_CLR_LIBS)
  # Setup references for /CLR
  set_target_properties(TestIPhreeqc PROPERTIES VS_DOTNET_REFERENCES "${CMAKE_CURRENT_SOURCE_DIR}/../src/phreeqcpp/ZedGraph.dll")
endif()

gtest_discover_tests(TestIPhreeqc)

add_executable(TestSelectedOutput TestSelectedOutput.cpp)
target_link_libraries(TestSelectedOutput IPhreeqc gtest gtest_main gmock)
gtest_discover_tests(TestSelectedOutput)

if (MSVC AND BUILD_SHARED_LIBS)
  # copy dlls for TestCVar
  add_custom_command(TARGET TestCVar PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:IPhreeqc> $<TARGET_FILE_DIR:TestCVar>
    )

  # copy dlls for TestVar
  add_custom_command(TARGET TestVar PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:IPhreeqc> $<TARGET_FILE_DIR:TestVar>
    )

  # copy dlls for TestIPhreeqcLib
  add_custom_command(TARGET TestIPhreeqcLib PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:IPhreeqc> $<TARGET_FILE_DIR:TestIPhreeqcLib>
    )

  # copy dlls for TestIPhreeqc
  add_custom_command(TARGET TestIPhreeqc PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:IPhreeqc> $<TARGET_FILE_DIR:TestIPhreeqc>
    )

  # copy dlls for TestSelectedOutput
  add_custom_command(TARGET TestSelectedOutput PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:IPhreeqc> $<TARGET_FILE_DIR:TestSelectedOutput>
    )
endif()


# databases
configure_file(
  phreeqc.dat.old
  phreeqc.dat.old
  COPYONLY
  )

configure_file(
  phreeqc.dat.90a6449
  phreeqc.dat.90a6449
  COPYONLY
  )

configure_file(
  ../database/phreeqc.dat
  phreeqc.dat
  COPYONLY
  )

configure_file(
  ../database/pitzer.dat
  pitzer.dat
  COPYONLY
  )

configure_file(
  ../database/wateq4f.dat
  wateq4f.dat
  COPYONLY
  )

configure_file(
  ../database/Amm.dat
  Amm.dat
  COPYONLY
  )

configure_file(
  ../database/frezchem.dat
  frezchem.dat
  COPYONLY
  )

configure_file(
  ../database/iso.dat
  iso.dat
  COPYONLY
  )

configure_file(
  ../database/llnl.dat
  llnl.dat
  COPYONLY
  )

configure_file(
  ../database/minteq.dat
  minteq.dat
  COPYONLY
  )

configure_file(
  ../database/minteq.v4.dat
  minteq.v4.dat
  COPYONLY
  )

configure_file(
  ../database/sit.dat
  sit.dat
  COPYONLY
  )

configure_file(
  ../database/ColdChem.dat
  ColdChem.dat
  COPYONLY
  )

configure_file(
  ../database/core10.dat
  core10.dat
  COPYONLY
  )

configure_file(
  ../database/Tipping_Hurley.dat
  Tipping_Hurley.dat
  COPYONLY
  )

configure_file(
  missing_e.dat
  missing_e.dat
  COPYONLY
  )

# misc files

configure_file(
  conv_fail.in
  conv_fail.in
  COPYONLY
  )

configure_file(
  dump
  dump
  COPYONLY
  )

configure_file(
  multi_punch
  multi_punch
  COPYONLY
  )

configure_file(
  multi_punch_no_set
  multi_punch_no_set
  COPYONLY
  )

configure_file(
  kinn20140218
  kinn20140218
  COPYONLY
  )

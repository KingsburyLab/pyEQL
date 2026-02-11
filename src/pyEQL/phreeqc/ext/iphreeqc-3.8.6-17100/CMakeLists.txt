# CMAKE_MSVC_RUNTIME_LIBRARY requires 3.15
cmake_minimum_required(VERSION 3.20)

project (
  IPhreeqc
  VERSION 3.8.6
  LANGUAGES CXX C
  )

# check if this is the root project
if (CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  set(STANDALONE_BUILD 1)
else()
  set(STANDALONE_BUILD 0)
endif()

if (STANDALONE_BUILD)
  # Set a default build type if none was specified
  set(default_build_type "Release")
  if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
      "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
  endif()
endif()

# override docdir on windows
if (WIN32 AND NOT CMAKE_INSTALL_DOCDIR)
  set(CMAKE_INSTALL_DOCDIR "" CACHE PATH "documentation root (doc)")
  set(CMAKE_INSTALL_DOCDIR "doc")
endif()

# set standard directory locations
include(GNUInstallDirs)

if (STANDALONE_BUILD AND MSVC)
  option (IPHREEQC_STATIC_RUNTIME "Build with a static runtime" OFF)
  if (IPHREEQC_STATIC_RUNTIME)
    # compile with static runtime
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
  endif()
endif()

# Fortran
option (IPHREEQC_ENABLE_MODULE "Enable Fortran module" ON)
if (STANDALONE_BUILD)
  if (IPHREEQC_FORTRAN_TESTING)
    enable_language(Fortran)
  else()
    if (CMAKE_Fortran_COMPILER)
      option (IPHREEQC_FORTRAN_TESTING "Build Fortran test" ON)
    else()
      option (IPHREEQC_FORTRAN_TESTING "Build Fortran test" OFF)
    endif()
  endif()
endif()

# compile Var.c as c++
set_source_files_properties(src/Var.c PROPERTIES LANGUAGE CXX)

if (STANDALONE_BUILD)
  option(BUILD_SHARED_LIBS "Build Shared Libraries" OFF)
  set(LIB_TYPE STATIC)
endif()

if (STANDALONE_BUILD)
  if (MSVC)
    option(BUILD_CLR_LIBS "Build CLR Libraries" OFF)
  endif()
endif()

if (BUILD_SHARED_LIBS)
  set(LIB_TYPE SHARED)
endif()

add_library(IPhreeqc ${LIB_TYPE} src/IPhreeqc.cpp)

target_sources(IPhreeqc
  PRIVATE
    src/CSelectedOutput.cpp
    src/CSelectedOutput.hxx
    src/CVar.hxx
    src/Debug.h
    src/ErrorReporter.hxx
    src/IPhreeqc.h
    src/IPhreeqc.hpp
    src/IPhreeqc_interface_F.cpp
    src/IPhreeqcCallbacks.h
    src/IPhreeqcLib.cpp
    src/phreeqcpp/advection.cpp
    src/phreeqcpp/basicsubs.cpp
    src/phreeqcpp/cl1.cpp
    src/phreeqcpp/common/Parser.cxx
    src/phreeqcpp/common/Parser.h
    src/phreeqcpp/common/PHRQ_base.cxx
    src/phreeqcpp/common/PHRQ_base.h
    src/phreeqcpp/common/PHRQ_exports.h
    src/phreeqcpp/common/PHRQ_io.cpp
    src/phreeqcpp/common/PHRQ_io.h
    src/phreeqcpp/common/phrqtype.h
    src/phreeqcpp/common/Utils.cxx
    src/phreeqcpp/common/Utils.h
    src/phreeqcpp/cvdense.cpp
    src/phreeqcpp/cvdense.h
    src/phreeqcpp/cvode.cpp
    src/phreeqcpp/cvode.h
    src/phreeqcpp/cxxKinetics.cxx
    src/phreeqcpp/cxxKinetics.h
    src/phreeqcpp/cxxMix.cxx
    src/phreeqcpp/cxxMix.h
    src/phreeqcpp/dense.cpp
    src/phreeqcpp/dense.h
    src/phreeqcpp/Dictionary.cpp
    src/phreeqcpp/Dictionary.h
    src/phreeqcpp/dumper.cpp
    src/phreeqcpp/dumper.h
    src/phreeqcpp/Exchange.cxx
    src/phreeqcpp/Exchange.h
    src/phreeqcpp/ExchComp.cxx
    src/phreeqcpp/ExchComp.h
    src/phreeqcpp/GasComp.cxx
    src/phreeqcpp/GasComp.h
    src/phreeqcpp/gases.cpp
    src/phreeqcpp/GasPhase.cxx
    src/phreeqcpp/GasPhase.h
    src/phreeqcpp/global_structures.h
    src/phreeqcpp/input.cpp
    src/phreeqcpp/integrate.cpp
    src/phreeqcpp/inverse.cpp
    src/phreeqcpp/ISolution.cxx
    src/phreeqcpp/ISolution.h
    src/phreeqcpp/ISolutionComp.cxx
    src/phreeqcpp/ISolutionComp.h
    src/phreeqcpp/isotopes.cpp
    src/phreeqcpp/kinetics.cpp
    src/phreeqcpp/KineticsComp.cxx
    src/phreeqcpp/KineticsComp.h
    src/phreeqcpp/mainsubs.cpp
    src/phreeqcpp/model.cpp
    src/phreeqcpp/NA.h
    src/phreeqcpp/NameDouble.cxx
    src/phreeqcpp/NameDouble.h
    src/phreeqcpp/NumKeyword.cxx
    src/phreeqcpp/NumKeyword.h
    src/phreeqcpp/nvector.cpp
    src/phreeqcpp/nvector.h
    src/phreeqcpp/nvector_serial.cpp
    src/phreeqcpp/nvector_serial.h
    src/phreeqcpp/parse.cpp
    src/phreeqcpp/PBasic.cpp
    src/phreeqcpp/PBasic.h
    src/phreeqcpp/phqalloc.cpp
    src/phreeqcpp/phqalloc.h
    src/phreeqcpp/Phreeqc.cpp
    src/phreeqcpp/Phreeqc.h
    src/phreeqcpp/PhreeqcKeywords/Keywords.cpp
    src/phreeqcpp/PhreeqcKeywords/Keywords.h
    src/phreeqcpp/PHRQ_io_output.cpp
    src/phreeqcpp/pitzer.cpp
    src/phreeqcpp/pitzer_structures.cpp
    src/phreeqcpp/PPassemblage.cxx
    src/phreeqcpp/PPassemblage.h
    src/phreeqcpp/PPassemblageComp.cxx
    src/phreeqcpp/PPassemblageComp.h
    src/phreeqcpp/prep.cpp
    src/phreeqcpp/Pressure.cxx
    src/phreeqcpp/Pressure.h
    src/phreeqcpp/print.cpp
    src/phreeqcpp/Reaction.cxx
    src/phreeqcpp/Reaction.h
    src/phreeqcpp/read.cpp
    src/phreeqcpp/ReadClass.cxx
    src/phreeqcpp/readtr.cpp
    src/phreeqcpp/runner.cpp
    src/phreeqcpp/runner.h
    src/phreeqcpp/SelectedOutput.cpp
    src/phreeqcpp/SelectedOutput.h
    src/phreeqcpp/Serializer.cxx
    src/phreeqcpp/Serializer.h
    src/phreeqcpp/sit.cpp
    src/phreeqcpp/smalldense.cpp
    src/phreeqcpp/smalldense.h
    src/phreeqcpp/Solution.cxx
    src/phreeqcpp/Solution.h
    src/phreeqcpp/SolutionIsotope.cxx
    src/phreeqcpp/SolutionIsotope.h
    src/phreeqcpp/spread.cpp
    src/phreeqcpp/SS.cxx
    src/phreeqcpp/SS.h
    src/phreeqcpp/SSassemblage.cxx
    src/phreeqcpp/SSassemblage.h
    src/phreeqcpp/SScomp.cxx
    src/phreeqcpp/SScomp.h
    src/phreeqcpp/step.cpp
    src/phreeqcpp/StorageBin.cxx
    src/phreeqcpp/StorageBin.h
    src/phreeqcpp/StorageBinList.cpp
    src/phreeqcpp/StorageBinList.h
    src/phreeqcpp/structures.cpp
    src/phreeqcpp/sundialsmath.cpp
    src/phreeqcpp/sundialsmath.h
    src/phreeqcpp/sundialstypes.h
    src/phreeqcpp/Surface.cxx
    src/phreeqcpp/Surface.h
    src/phreeqcpp/SurfaceCharge.cxx
    src/phreeqcpp/SurfaceCharge.h
    src/phreeqcpp/SurfaceComp.cxx
    src/phreeqcpp/SurfaceComp.h
    src/phreeqcpp/System.cxx
    src/phreeqcpp/System.h
    src/phreeqcpp/tally.cpp
    src/phreeqcpp/Temperature.cxx
    src/phreeqcpp/Temperature.h
    src/phreeqcpp/tidy.cpp
    src/phreeqcpp/transport.cpp
    src/phreeqcpp/Use.cpp
    src/phreeqcpp/Use.h
    src/phreeqcpp/UserPunch.cpp
    src/phreeqcpp/UserPunch.h
    src/phreeqcpp/utilities.cpp
    src/thread.h
    src/Var.c
    src/Var.h
    src/Version.h
)

target_include_directories(IPhreeqc
  PUBLIC
    $<BUILD_INTERFACE:${IPhreeqc_SOURCE_DIR}/src>
    $<BUILD_INTERFACE:${IPhreeqc_SOURCE_DIR}/src/phreeqcpp>
    $<BUILD_INTERFACE:${IPhreeqc_SOURCE_DIR}/src/phreeqcpp/common>
    $<BUILD_INTERFACE:${IPhreeqc_SOURCE_DIR}/src/phreeqcpp/PhreeqcKeywords>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  )

target_compile_definitions(IPhreeqc PRIVATE SWIG_SHARED_OBJ)
target_compile_definitions(IPhreeqc PRIVATE USE_PHRQ_ALLOC)

if (NOT IPHREEQC_ENABLE_MODULE)
  target_compile_definitions(IPhreeqc
    PUBLIC
      IPHREEQC_NO_FORTRAN_MODULE
  )

  target_sources(IPhreeqc
    PRIVATE
      src/fimpl.h
      src/fwrap.cpp
      src/fwrap1.cpp
      src/fwrap2.cpp
      src/fwrap3.cpp
      src/fwrap4.cpp
      src/fwrap5.cpp
      src/fwrap6.cpp
      src/fwrap7.cpp
      src/fwrap8.cpp
  )
endif()

if (${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
  target_compile_options(IPhreeqc PRIVATE /wd4251 /wd4275 /wd4793)
endif()

# Disable deprecation warnings for standard C and STL functions in VS2005
# and later
if (MSVC_VERSION EQUAL 1400 OR MSVC_VERSION GREATER 1400)
  target_compile_definitions(IPhreeqc PUBLIC _CRT_SECURE_NO_DEPRECATE)
  target_compile_definitions(IPhreeqc PUBLIC _SCL_SECURE_NO_WARNINGS)
endif()

if (WIN32 AND BUILD_SHARED_LIBS)
  # Required to build IMPLIB
  # (Seems to be automatically set when using Visual Studio as the generator)
  target_compile_definitions(IPhreeqc PRIVATE _WINDLL)
endif()

# Allow user to override POSTFIX properties (but mandate them so that
# all builds can be installed to the same directory)
if (NOT CMAKE_DEBUG_POSTFIX)
  set_target_properties(IPhreeqc PROPERTIES DEBUG_POSTFIX "d")
endif()

if (NOT CMAKE_MINSIZEREL_POSTFIX)
  set_target_properties(IPhreeqc PROPERTIES MINSIZEREL_POSTFIX "msr")
endif()

if (NOT CMAKE_RELWITHDEBINFO_POSTFIX)
  set_target_properties(IPhreeqc PROPERTIES RELWITHDEBINFO_POSTFIX "rwd")
endif()

set(IPhreeqc_Headers
  ${PROJECT_SOURCE_DIR}/src/IPhreeqc.h
  ${PROJECT_SOURCE_DIR}/src/IPhreeqc.hpp
  ${PROJECT_SOURCE_DIR}/src/IPhreeqcCallbacks.h
  ${PROJECT_SOURCE_DIR}/src/phreeqcpp/PhreeqcKeywords/Keywords.h
  ${PROJECT_SOURCE_DIR}/src/phreeqcpp/common/PHRQ_exports.h
  ${PROJECT_SOURCE_DIR}/src/phreeqcpp/common/PHRQ_io.h
  ${PROJECT_SOURCE_DIR}/src/Var.h
  )

if (IPHREEQC_ENABLE_MODULE)
  set(IPhreeqc_SRC
    ${PROJECT_SOURCE_DIR}/src/IPhreeqc_interface.F90
    ${PROJECT_SOURCE_DIR}/src/README.Fortran
    )
else()
  # always install README.Fortran
  set(IPhreeqc_SRC
    ${PROJECT_SOURCE_DIR}/src/README.Fortran
    )
  # install old fortran include files
  set(IPhreeqc_Headers
    ${IPhreeqc_Headers}
    ${PROJECT_SOURCE_DIR}/src/IPhreeqc.f.inc
    ${PROJECT_SOURCE_DIR}/src/IPhreeqc.f90.inc
    )
endif()

# Setup references for /CLR
if (MSVC AND BUILD_SHARED_LIBS AND BUILD_CLR_LIBS)
  if (MSVC_VERSION LESS 1600)
    message(FATAL_ERROR "CLR options must be set manually for versions prior to Visual Studio 2010")
  endif()
  set_target_properties(IPhreeqc PROPERTIES VS_DOTNET_REFERENCES "System;System.Drawing;System.Windows.Forms;${CMAKE_CURRENT_SOURCE_DIR}/src/phreeqcpp/ZedGraph.dll")
  set_target_properties(IPhreeqc PROPERTIES COMMON_LANGUAGE_RUNTIME "")
  target_compile_definitions(IPhreeqc PRIVATE "MULTICHART")
endif()

set_target_properties(IPhreeqc PROPERTIES
  VERSION "${IPhreeqc_VERSION}"
  SOVERSION "${IPhreeqc_VERSION_MAJOR}"
  )

# install
if (STANDALONE_BUILD)

  include(CMakePackageConfigHelpers)

  configure_package_config_file(IPhreeqcConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/IPhreeqcConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/IPhreeqc
    )

  write_basic_package_version_file("IPhreeqcConfigVersion.cmake"
    VERSION ${IPhreeqc_VERSION}
    COMPATIBILITY SameMajorVersion
    )

  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/IPhreeqcConfig.cmake ${CMAKE_CURRENT_BINARY_DIR}/IPhreeqcConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/IPhreeqc
    )

  install(TARGETS  IPhreeqc
    EXPORT   IPhreeqcTargets
    LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    )
  install(EXPORT IPhreeqcTargets
    FILE IPhreeqcTargets.cmake
    NAMESPACE IPhreeqc::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/IPhreeqc
    )

  install(FILES ${IPhreeqc_Headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  if (WIN32)
    install(FILES ${IPhreeqc_SRC} DESTINATION src)
  else()
    install(FILES ${IPhreeqc_SRC} DESTINATION ${CMAKE_INSTALL_DOCDIR}/src)
  endif()
  if (BUILD_CLR_LIBS)
    install(FILES "${PROJECT_SOURCE_DIR}/src/phreeqcpp/ZedGraph.dll" DESTINATION ${CMAKE_INSTALL_BINDIR})
  endif()
endif()

if (BUILD_CLR_LIBS)
  if (NOT BUILD_SHARED_LIBS)
    message(FATAL_ERROR "BUILD_CLR_LIBS requires that BUILD_SHARED_LIBS be set")
  endif()

  if (IPHREEQC_STATIC_RUNTIME)
    message(FATAL_ERROR "BUILD_CLR_LIBS and IPHREEQC_STATIC_RUNTIME are mutually exclusive")
  endif()

  # CLR files
  target_sources(IPhreeqc
    PRIVATE
      src/phreeqcpp/ChartHandler.cpp
      src/phreeqcpp/ChartHandler.h
      src/phreeqcpp/ChartObject.cpp
      src/phreeqcpp/ChartObject.h
      src/phreeqcpp/CurveObject.cpp
      src/phreeqcpp/CurveObject.h
  )

endif()

include (CTest)

if (STANDALONE_BUILD)

  add_subdirectory(database)
  add_subdirectory(doc)
  add_subdirectory(examples)

  if (BUILD_TESTING)
    add_subdirectory(tests)
  endif()


  if (BUILD_TESTING)  # may need to add MSVC version check
    include(FetchContent)

    # Avoid warning about DOWNLOAD_EXTRACT_TIMESTAMP in CMake 3.24:
    if (CMAKE_VERSION VERSION_GREATER_EQUAL "3.24.0")
      cmake_policy(SET CMP0135 NEW)
    endif()    

    FetchContent_Declare(
      googletest
      URL https://github.com/google/googletest/archive/release-1.12.1.tar.gz
      )

    mark_as_advanced(
      FETCHCONTENT_BASE_DIR
      FETCHCONTENT_FULLY_DISCONNECTED
      FETCHCONTENT_QUIET
      FETCHCONTENT_SOURCE_DIR_GOOGLETEST
      FETCHCONTENT_UPDATES_DISCONNECTED
      FETCHCONTENT_UPDATES_DISCONNECTED_GOOGLETEST
      )

    # Prevent GoogleTest from overriding our compiler/linker options
    # when building with Visual Studio
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

    set(BUILD_GTEST ON CACHE BOOL "" FORCE)
    set(BUILD_GMOCK ON CACHE BOOL "" FORCE)

    set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
    set(INSTALL_GMOCK OFF CACHE BOOL "" FORCE)

    mark_as_advanced(
      BUILD_GTEST
      BUILD_GMOCK
      INSTALL_GTEST
      INSTALL_GMOCK
      gmock_build_tests
      gtest_build_tests
      gtest_build_samples
      gtest_disable_pthreads
      gtest_force_shared_crt
      gtest_hide_internal_symbols
      )

    FetchContent_GetProperties(googletest)
    if (NOT googletest_POPULATED)
      # Always build googletest static
      set(SAVE_BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})
      set(BUILD_SHARED_LIBS OFF)
      FetchContent_MakeAvailable(googletest)
      set(BUILD_SHARED_LIBS ${SAVE_BUILD_SHARED_LIBS})
    endif()

    add_subdirectory(gtest)
  endif()
endif()

# get_cmake_property(_variableNames VARIABLES)
# list (SORT _variableNames)
# foreach (_variableName ${_variableNames})
#   message(STATUS "${_variableName}=${${_variableName}}")
# endforeach()

cmake_minimum_required(VERSION 3.15...3.26)
project(
  "${SKBUILD_PROJECT_NAME}"
  LANGUAGES CXX
  VERSION "${SKBUILD_PROJECT_VERSION}")

option(PYEQL_ENABLE_EXT "Build and install IPhreeqc support" ON)

if(PYEQL_ENABLE_EXT)
    find_package(
      Python
      COMPONENTS Interpreter Development.Module
      REQUIRED)
    find_package(pybind11 CONFIG REQUIRED)

    # -------------------------------------------------------
    # Setup for IPhreeqc
    # -------------------------------------------------------
    set(IPHREEQC_BASE ${CMAKE_CURRENT_SOURCE_DIR}/src/pyEQL/phreeqc/ext/iphreeqc-3.8.6-17100)
    set(IPHREEQC_DIR         ${IPHREEQC_BASE})
    set(IPHREEQC_INCLUDE_DIRS
        ${IPHREEQC_BASE}/src
        ${IPHREEQC_BASE}/src/phreeqcpp/common
    )
    include_directories(${IPHREEQC_INCLUDE_DIRS})
    set(IPHREEQC_BUILD_DIR   ${IPHREEQC_BASE}/build)
    add_subdirectory(${IPHREEQC_DIR} ${IPHREEQC_BUILD_DIR})

    # Hide IPhreeqc's symbols when it is statically linked into _iphreeqc.so.
    # Upstream IPhreeqc does not build with -fvisibility=hidden, so without
    # this every Phreeqc/IPhreeqcLib symbol gets re-exported from _iphreeqc.so.
    # That collides with phreeqpython's own copy of IPhreeqc (a different
    # build), causing segfault once both modules are loaded in the same Python
    # process.
    # See https://linuxvox.com/blog/c-fvisibility-hidden-fvisibility-inlines-hidden/
    set_target_properties(IPhreeqc PROPERTIES
        POSITION_INDEPENDENT_CODE ON
        C_VISIBILITY_PRESET hidden
        CXX_VISIBILITY_PRESET hidden
        VISIBILITY_INLINES_HIDDEN ON)
    # -------------------------------------------------------

    # -------------------------------------------------------
    # Setup for pyEQL.phreeqc._iphreeqc
    # -------------------------------------------------------
    pybind11_add_module(_iphreeqc ${CMAKE_CURRENT_SOURCE_DIR}/src/pyEQL/phreeqc/bindings.cpp)
    target_compile_definitions(_iphreeqc
                               PRIVATE VERSION_INFO=${PROJECT_VERSION})
    target_link_libraries(_iphreeqc PRIVATE IPhreeqc)
    # std::variant in iphreeqc_wrapper.hpp (MSVC defaults to C++14)
    target_compile_features(_iphreeqc PRIVATE cxx_std_17)
    install(TARGETS _iphreeqc LIBRARY DESTINATION pyEQL/phreeqc COMPONENT python)

    # -------------------------------------------------------
    # IPhreeqc databases
    # -------------------------------------------------------
    # Till we find a way to extract artifacts from IPhreeqc's install folder
    # through scikit-build-core, we define a COMPONENT that allows us to
    # copy the database files.
    set(IPHREEQC_DATABASE_DIR   ${IPHREEQC_BASE}/database)
    install(DIRECTORY
           "${IPHREEQC_DATABASE_DIR}"
           DESTINATION pyEQL/phreeqc COMPONENT iphreeqc_database)
    # Make it a python module so we can access it using importlib.resources
    file(GENERATE
        OUTPUT "__init__.py"
        CONTENT ""
    )
    install(FILES
            "${CMAKE_CURRENT_BINARY_DIR}/__init__.py"
            DESTINATION pyEQL/phreeqc/database
            COMPONENT iphreeqc_database
    )
    # -------------------------------------------------------
endif()

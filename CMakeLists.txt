cmake_minimum_required(VERSION 3.15...3.26)
project(
  "${SKBUILD_PROJECT_NAME}"
  LANGUAGES CXX
  VERSION "${SKBUILD_PROJECT_VERSION}")

find_package(
  Python
  COMPONENTS Interpreter Development.Module
  REQUIRED)
find_package(pybind11 CONFIG REQUIRED)

# -------------------------------------------------------
# IPhreeqc setup
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
set_target_properties(IPhreeqc PROPERTIES POSITION_INDEPENDENT_CODE ON)
# -------------------------------------------------------

# -------------------------------------------------------
# Setup for pyEQL._phreeqc
# -------------------------------------------------------
pybind11_add_module(_phreeqc ${CMAKE_CURRENT_SOURCE_DIR}/src/pyEQL/phreeqc/bindings.cpp)
target_compile_definitions(_phreeqc
                           PRIVATE VERSION_INFO=${PROJECT_VERSION})
target_link_libraries(_phreeqc PRIVATE IPhreeqc)
install(TARGETS _phreeqc LIBRARY DESTINATION pyEQL)

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
# -------------------------------------------------------

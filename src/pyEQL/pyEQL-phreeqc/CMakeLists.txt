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
set(IPHREEQC_BASE ${CMAKE_CURRENT_SOURCE_DIR}/ext/iphreeqc-3.8.6-17100)
set(IPHREEQC_DIR         ${IPHREEQC_BASE})
set(IPHREEQC_INCLUDE_DIRS
    ${IPHREEQC_BASE}/src
    ${IPHREEQC_BASE}/src/phreeqcpp/common
)
set(IPHREEQC_BUILD_DIR   ${IPHREEQC_BASE}/build)
add_subdirectory(${IPHREEQC_DIR} ${IPHREEQC_BUILD_DIR})
set_target_properties(IPhreeqc PROPERTIES POSITION_INDEPENDENT_CODE ON)
# -------------------------------------------------------

# -------------------------------------------------------
# IPhreeqc databases
# -------------------------------------------------------
# Till we find a way to extract artifacts from IPhreeqc's install folder
# through scikit-build-core, we define a COMPONENT that allows us to
# copy the database files.
set(IPHREEQC_DATABASE_DIR   ${IPHREEQC_BASE}/database)
install(DIRECTORY
      "${IPHREEQC_DATABASE_DIR}"
      DESTINATION pyEQL_phreeqc COMPONENT iphreeqc_database)
# -------------------------------------------------------

include_directories(${IPHREEQC_INCLUDE_DIRS})

pybind11_add_module(_bindings ${CMAKE_CURRENT_SOURCE_DIR}/bindings.cpp)

target_compile_definitions(_bindings
                           PRIVATE VERSION_INFO=${PROJECT_VERSION})

target_link_libraries(_bindings PRIVATE IPhreeqc)

install(TARGETS _bindings LIBRARY DESTINATION pyEQL_phreeqc)

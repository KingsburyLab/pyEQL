EXTRA_DIST =\
	CMakeLists.txt\
	c/CMakeLists.txt\
	c/advect/CMakeLists.txt\
	c/advect/CMakeLists.txt.in\
	c/advect/README.txt\
	com/CMakeLists.txt\
	com/README.txt\
	com/excel/CMakeLists.txt\
	com/python/CMakeLists.txt\
	cpp/CMakeLists.txt\
	cpp/advect/CMakeLists.txt\
	cpp/advect/CMakeLists.txt.in\
	cpp/advect/README.txt\
	fortran/CMakeLists.txt\
	fortran/advect/CMakeLists.txt\
	fortran/advect/CMakeLists.txt.in\
	fortran/advect/README.txt\
	using-cmake/CMakeLists.txt\
	using-cmake/CMakeLists.txt.in\
	using-cmake/ex2\
	using-cmake/main.cpp\
	using-cmake/phreeqc.dat\
	using-cmake/post-install.cmake.in\
	using-cmake/README.txt

EXAMPLES_DIR = $(docdir)/examples

# c
#
example_cdir = $(EXAMPLES_DIR)/c
example_c_advectdir = $(EXAMPLES_DIR)/c/advect

dist_example_c_advect_DATA = $(c_advect)

c_advect= \
    c/advect/advect.c \
    c/advect/ic \
    c/advect/phreeqc.dat

# c++
#
example_cppdir = $(EXAMPLES_DIR)/cpp
example_cpp_advectdir = $(EXAMPLES_DIR)/cpp/advect

dist_example_cpp_advect_DATA = $(cpp_advect)

cpp_advect= \
    cpp/advect/advect.cpp \
    cpp/advect/ic \
    cpp/advect/phreeqc.dat

# fortran
#
example_fortrandir = $(EXAMPLES_DIR)/fortran
example_fortran_advectdir = $(EXAMPLES_DIR)/fortran/advect

dist_example_fortran_advect_DATA = $(fortran_advect)

fortran_advect= \
    fortran/advect/advect.F90 \
    fortran/advect/ic \
    fortran/advect/phreeqc.dat

# com excel
#
example_comdir = $(EXAMPLES_DIR)/com
example_com_exceldir = $(EXAMPLES_DIR)/com/excel

dist_example_com_excel_DATA = $(com_excel)

com_excel= \
    com/excel/phreeqc.dat \
    com/excel/runphreeqc.xls \
    com/excel/withcallback.xls

# com python
#
example_com_pythondir = $(EXAMPLES_DIR)/com/python

dist_example_com_python_DATA = $(com_python)

com_python= \
    com/python/Gypsum.py \
    com/python/parallel_advect.py \
    com/python/phreeqc.dat \
    com/python/pitzer.dat \
    com/python/wateq4f.dat

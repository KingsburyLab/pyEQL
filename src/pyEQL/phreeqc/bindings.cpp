#include <pybind11/pybind11.h>
#include <pybind11/stl.h>  // std::variant/std::string conversions
#include "iphreeqc_wrapper.hpp"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

PYBIND11_MODULE(_iphreeqc, m) {

    py::class_<IPhreeqcWrapper>(m, "IPhreeqcWrapper")
        .def(py::init<>())
        .def("load_database", &IPhreeqcWrapper::load_database)
        .def("run_string", &IPhreeqcWrapper::run_string)
        .def("get_error_string", &IPhreeqcWrapper::get_error_string)
        .def("get_selected_output_row_count", &IPhreeqcWrapper::get_selected_output_row_count)
        .def("get_selected_output_column_count", &IPhreeqcWrapper::get_selected_output_column_count)
        .def("get_value", &IPhreeqcWrapper::get_value)
        .def("get_component_count", &IPhreeqcWrapper::get_component_count)
        .def("get_component", &IPhreeqcWrapper::get_component)
        .def("set_dump_string_on", &IPhreeqcWrapper::set_dump_string_on)
        .def("get_dump_string", &IPhreeqcWrapper::get_dump_string)
        .def("set_log_string_on", &IPhreeqcWrapper::set_log_string_on)
        .def("get_log_string", &IPhreeqcWrapper::get_log_string);


#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}

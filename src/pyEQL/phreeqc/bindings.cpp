#include <pybind11/pybind11.h>
#include "iphreeqc_wrapper.cpp"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

// managed Var that doesn't leak memory
class PyVar {
    public: VAR var;
public:
    PyVar() {
        VarInit(&this->var);
    }

    ~PyVar() {
        VarClear(&this->var);
    }
};


PYBIND11_MODULE(_phreeqc, m) {

    py::enum_<VAR_TYPE>(m, "PY_VAR_TYPE", py::module_local())
        .value("TT_EMPTY", TT_EMPTY)
        .value("TT_ERROR", TT_ERROR)
        .value("TT_LONG", TT_LONG)
        .value("TT_DOUBLE", TT_DOUBLE)
        .value("TT_STRING", TT_STRING);

    py::enum_<VRESULT>(m, "PY_VRESULT", py::module_local())
        .value("VR_OK", VR_OK)
        .value("VR_OUTOFMEMORY", VR_OUTOFMEMORY)
        .value("VR_BADVARTYPE", VR_BADVARTYPE)
        .value("VR_INVALIDARG", VR_INVALIDARG)
        .value("VR_INVALIDROW", VR_INVALIDROW)
        .value("VR_INVALIDCOL", VR_INVALIDCOL);

    // The VAR type is exposed to Python, but only indirectly through an
    // attribute of the PyVar class, which is why we name the class
    // "_VAR".
    py::class_<VAR>(m, "_VAR")
        .def_readwrite("type", &VAR::type)
        .def_readwrite("lVal", &VAR::lVal)
        .def_readwrite("dVal", &VAR::dVal)
        .def_property("sVal",   // property since str needs to be on the heap
            [](const VAR& v) {  // getter
                return std::string(v.sVal ? v.sVal : "");
            },
            [](VAR& v, const std::string& val) {  // setter
                v.sVal = strdup(val.c_str());  // needs to be freed, which is done in VarClear
            })
        .def_readwrite("vresult", &VAR::vresult);

    py::class_<PyVar>(m, "PyVar")
        .def(py::init<>())
        // For the following readwrite attribute to work, pybind11 needs to
        // know how to expose VAR as a Python object, which is why we defined
        // the class above.
        .def_readwrite("var", &PyVar::var);

    py::class_<IPhreeqcWrapper>(m, "PyIPhreeqc")
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

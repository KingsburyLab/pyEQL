#pragma once

#include <IPhreeqc.h>
#include <string>
#include <stdexcept>
#include <variant>

// Value of one selected-output cell; pybind11 converts this to None/int/float/str.
using CellValue = std::variant<std::monostate, long, double, std::string>;

class IPhreeqcWrapper {
public:
    IPhreeqcWrapper() {
        id = CreateIPhreeqc();
        if (id < 0) throw std::runtime_error("Failed to create IPhreeqc instance");
    }

    ~IPhreeqcWrapper() {
        DestroyIPhreeqc(id);
    }

    void load_database(const std::string& db_path) {
        if (LoadDatabase(id, db_path.c_str()) != 0) {
            throw std::runtime_error(get_error_string());
        }
    }

    void run_string(const std::string& input) {
        if (RunString(id, input.c_str()) != 0) {
            throw std::runtime_error(get_error_string());
        }
    }

    std::string get_error_string() {
        const char* err = GetErrorString(id);
        return std::string(err ? err : "");
    }

    int get_selected_output_row_count() {
        return GetSelectedOutputRowCount(id);
    }

    int get_selected_output_column_count() {
        return GetSelectedOutputColumnCount(id);
    }

    CellValue get_value(int row, int col) {
        // VAR is a C tagged union whose sVal is heap-allocated, so it must be
        // paired with VarClear. Decode it here so it never reaches Python.
        VAR var;
        VarInit(&var);
        if (GetSelectedOutputValue(id, row, col, &var) != VR_OK) {
            VarClear(&var);
            throw std::runtime_error(get_error_string());
        }
        CellValue value;
        switch (var.type) {
            case TT_LONG:   value = var.lVal; break;
            case TT_DOUBLE: value = var.dVal; break;
            case TT_STRING: value = std::string(var.sVal ? var.sVal : ""); break;
            default:        break;  // TT_EMPTY
        }
        VarClear(&var);
        return value;
    }

    int get_component_count() {
        return GetComponentCount(id);
    }

    std::string get_component(int i) {
        const char* s = GetComponent(id, i);
        return std::string(s ? s : "");
    }

    std::string get_dump_string() {
        const char* s = GetDumpString(id);
        return std::string(s ? s : "");
    }

    void set_dump_string_on(int i) {
        if (SetDumpStringOn(id, i) != 0) {
            throw std::runtime_error(get_error_string());
        }
    }

    std::string get_log_string() {
        const char* s = GetLogString(id);
        return std::string(s ? s : "");
    }

    void set_log_string_on(int i) {
        if (SetLogStringOn(id, i) != 0) {
            throw std::runtime_error(get_error_string());
        }
    }

private:
    int id;
};

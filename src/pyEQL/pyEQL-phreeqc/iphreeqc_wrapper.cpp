#include <IPhreeqc.h>
#include <string>
#include <vector>
#include <stdexcept>

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

    int get_value(int row, int col, VAR& var) {
        if (GetSelectedOutputValue(id, row, col, &var) != VR_OK) {
            throw std::runtime_error(get_error_string());
        }
        return var.vresult;
    }

    int get_component_count() {
        return GetComponentCount(id);
    }

    std::string get_component(int i) {
        return GetComponent(id, i);
    }

    std::string get_dump_string() {
        return GetDumpString(id);
    }

    int set_dump_string_on(int i) {
        return SetDumpStringOn(id, i);
    }

    std::string get_log_string() {
        return GetLogString(id);
    }

    int set_log_string_on(int i) {
        return SetLogStringOn(id, i);
    }

private:
    int id;
};

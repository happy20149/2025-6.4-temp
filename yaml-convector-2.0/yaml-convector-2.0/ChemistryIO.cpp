// ChemistryIO.cpp
#include "ChemistryIO.h"
#include <iostream>
#include <fstream>
#include <stdexcept>

#ifdef HAVE_YAML_CPP
// YAML::Node to YamlValue conversion
ChemistryIO::YamlValue::YamlValue(const YAML::Node& node) {
    if (node.IsNull()) {
        m_type = Type::Null;
    }
    else if (node.IsScalar()) {
        try {
            if (node.as<std::string>() == "true" || node.as<std::string>() == "false") {
                m_type = Type::Boolean;
                m_bool = node.as<bool>();
            }
            else {
                try {
                    double number = node.as<double>();
                    m_type = Type::Number;
                    m_number = number;
                }
                catch (...) {
                    m_type = Type::String;
                    m_string = node.as<std::string>();
                }
            }
        }
        catch (...) {
            m_type = Type::String;
            m_string = node.as<std::string>();
        }
    }
    else if (node.IsMap()) {
        m_type = Type::Map;
        for (const auto& kv : node) {
            m_map[kv.first.as<std::string>()] = YamlValue(kv.second);
        }
    }
    else if (node.IsSequence()) {
        m_type = Type::Sequence;
        for (const auto& item : node) {
            m_sequence.push_back(YamlValue(item));
        }
    }
}
#endif

std::string ChemistryIO::YamlValue::asString() const {
    if (m_type != Type::String) {
        throw std::runtime_error("Value is not a string");
    }
    return m_string;
}

double ChemistryIO::YamlValue::asNumber() const {
    if (m_type != Type::Number) {
        throw std::runtime_error("Value is not a number");
    }
    return m_number;
}

bool ChemistryIO::YamlValue::asBoolean() const {
    if (m_type != Type::Boolean) {
        throw std::runtime_error("Value is not a boolean");
    }
    return m_bool;
}

const std::map<std::string, ChemistryIO::YamlValue>& ChemistryIO::YamlValue::asMap() const {
    if (m_type != Type::Map) {
        throw std::runtime_error("Value is not a map");
    }
    return m_map;
}

const std::vector<ChemistryIO::YamlValue>& ChemistryIO::YamlValue::asSequence() const {
    if (m_type != Type::Sequence) {
        throw std::runtime_error("Value is not a sequence");
    }
    return m_sequence;
}

void ChemistryIO::YamlValue::print(int indent) const {
    std::string indentStr(indent, ' ');

    switch (m_type) {
    case Type::Null:
        std::cout << indentStr << "null" << std::endl;
        break;
    case Type::String:
        std::cout << indentStr << "\"" << m_string << "\"" << std::endl;
        break;
    case Type::Number:
        std::cout << indentStr << m_number << std::endl;
        break;
    case Type::Boolean:
        std::cout << indentStr << (m_bool ? "true" : "false") << std::endl;
        break;
    case Type::Map:
        std::cout << indentStr << "{" << std::endl;
        for (const auto& pair : m_map) {
            const auto& key = pair.first;
            const auto& value = pair.second;
            std::cout << indentStr << "  " << key << ": ";
            value.print(indent + 4);
        }
        std::cout << indentStr << "}" << std::endl;
        break;
    case Type::Sequence:
        std::cout << indentStr << "[" << std::endl;
        for (const auto& item : m_sequence) {
            std::cout << indentStr << "  - ";
            item.print(indent + 4);
        }
        std::cout << indentStr << "]" << std::endl;
        break;
    }
}

// File loading functions
ChemistryIO::YamlValue ChemistryIO::loadFile(const std::string& filename) {
#ifdef HAVE_YAML_CPP
    try {
        YAML::Node node = YAML::LoadFile(filename);
        return YamlValue(node);
    }
    catch (const std::exception& e) {
        std::cerr << "Failed to load YAML file: " << e.what() << std::endl;
        throw;
    }
#else
    std::cerr << "YAML-cpp support not available. Cannot load file: " << filename << std::endl;
    return YamlValue(); // Return null value
#endif
}

ChemistryIO::YamlValue ChemistryIO::loadString(const std::string& yaml) {
#ifdef HAVE_YAML_CPP
    try {
        YAML::Node node = YAML::Load(yaml);
        return YamlValue(node);
    }
    catch (const std::exception& e) {
        std::cerr << "Failed to parse YAML string: " << e.what() << std::endl;
        throw;
    }
#else
    std::cerr << "YAML-cpp support not available. Cannot parse YAML string." << std::endl;
    return YamlValue(); // Return null value
#endif
}

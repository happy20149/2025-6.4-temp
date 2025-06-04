// ChemistryIO.cpp
#include "ChemistryIO.h"
#include <iostream>
#include <fstream>
#include <stdexcept>

// 从YamlValue.cpp移植的代码

// 将YAML::Node转换为YamlValue
ChemistryIO::YamlValue::YamlValue(const YAML::Node& node) {
    if (node.IsNull()) {
        m_type = Type::Null;
    }
    else if (node.IsScalar()) {
        // 尝试转换为不同的标量类型
        try {
            // 尝试转换为布尔值
            if (node.as<std::string>() == "true" || node.as<std::string>() == "false") {
                m_type = Type::Boolean;
                m_bool = node.as<bool>();
            }
            // 尝试转换为数字
            else {
                try {
                    double number = node.as<double>();
                    m_type = Type::Number;
                    m_number = number;
                }
                catch (...) {
                    // 如果转换失败，默认为字符串
                    m_type = Type::String;
                    m_string = node.as<std::string>();
                }
            }
        }
        catch (...) {
            // 如果转换失败，默认为字符串
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

std::string ChemistryIO::YamlValue::asString() const {
    if (m_type != Type::String) {
        throw std::runtime_error("尝试将非字符串类型作为字符串访问");
    }
    return m_string;
}

double ChemistryIO::YamlValue::asNumber() const {
    if (m_type != Type::Number) {
        throw std::runtime_error("尝试将非数字类型作为数字访问");
    }
    return m_number;
}

bool ChemistryIO::YamlValue::asBoolean() const {
    if (m_type != Type::Boolean) {
        throw std::runtime_error("尝试将非布尔类型作为布尔值访问");
    }
    return m_bool;
}

const std::map<std::string, ChemistryIO::YamlValue>& ChemistryIO::YamlValue::asMap() const {
    if (m_type != Type::Map) {
        throw std::runtime_error("尝试将非映射表类型作为映射表访问");
    }
    return m_map;
}

const std::vector<ChemistryIO::YamlValue>& ChemistryIO::YamlValue::asSequence() const {
    if (m_type != Type::Sequence) {
        throw std::runtime_error("尝试将非序列类型作为序列访问");
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
        //for (const auto& [key, value] : m_map) {
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

// 从YamlParser.cpp移植的文件加载功能
ChemistryIO::YamlValue ChemistryIO::loadFile(const std::string& filename) {
    try {
        YAML::Node node = YAML::LoadFile(filename);
        return YamlValue(node);
    }
    catch (const std::exception& e) {
        std::cerr << "加载YAML文件失败: " << e.what() << std::endl;
        throw;
    }
}

ChemistryIO::YamlValue ChemistryIO::loadString(const std::string& yaml) {
    try {
        YAML::Node node = YAML::Load(yaml);
        return YamlValue(node);
    }
    catch (const std::exception& e) {
        std::cerr << "解析YAML字符串失败: " << e.what() << std::endl;
        throw;
    }
}


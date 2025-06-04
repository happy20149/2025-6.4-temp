// ChemistryIO.h
#pragma once
#include <string>
#include <map>
#include <vector>
//#include <any>

#ifdef HAVE_YAML_CPP
#include <yaml-cpp/yaml.h>
#endif

// 从YamlParser.h移植的类和函数
class ChemistryIO {
public:
    // YamlValue类定义
    class YamlValue {
    public:
        //map:映射表(键值对集合)  sequence:序列(元素集合)
        enum class Type {
            Null, String, Number, Boolean, Map, Sequence
        };        // 构造函数
        YamlValue() : m_type(Type::Null) {}
        YamlValue(const std::string& value) : m_type(Type::String), m_string(value) {}
        YamlValue(double value) : m_type(Type::Number), m_number(value) {}
        YamlValue(bool value) : m_type(Type::Boolean), m_bool(value) {}
#ifdef HAVE_YAML_CPP
        YamlValue(const YAML::Node& node);
#endif

        // 类型检查
        bool isNull() const { return m_type == Type::Null; }
        bool isString() const { return m_type == Type::String; }
        bool isNumber() const { return m_type == Type::Number; }
        bool isBoolean() const { return m_type == Type::Boolean; }
        bool isMap() const { return m_type == Type::Map; }
        bool isSequence() const { return m_type == Type::Sequence; }

        // 值访问
        std::string asString() const;
        double asNumber() const;
        bool asBoolean() const;
        const std::map<std::string, YamlValue>& asMap() const;
        const std::vector<YamlValue>& asSequence() const;

        // 打印值
        void print(int indent = 0) const;

    private:
        Type m_type;
        std::string m_string;
        double m_number = 0.0;
        bool m_bool = false;
        std::map<std::string, YamlValue> m_map;
        std::vector<YamlValue> m_sequence;
    };

    // YamlParser静态函数
    static YamlValue loadFile(const std::string& filename);
    static YamlValue loadString(const std::string& yaml);
};
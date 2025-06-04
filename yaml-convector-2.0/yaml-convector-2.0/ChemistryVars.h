
#pragma once
#include <string>
#include <vector>
#include <map>


class ChemistryVars {
public:
    // 反应数据结构
    struct ReactionData {
        std::string equation;
        std::string type;

        // 速率常数
        struct {
            double A = 0.0;
            std::string A_units;
            double b = 0.0;
            double Ea = 0.0;
            std::string Ea_units;

            // PLOG支持
            struct PLOGPoint {
                double pressure;  // atm
                double A;
                double b;
                double Ea;
            };

            bool isPressureDependent = false;
            std::vector<PLOGPoint> plogData;
        } rateConstant;

        // 第三体效率
        std::map<std::string, double> efficiencies;

        // 低压极限
        struct {
            double A = 0.0;
            double b = 0.0;
            double Ea = 0.0;
        } lowPressure;

        // Troe参数
        struct {
            double a = 0.0;
            double T_star = 0.0;
            double T_double_star = 0.0;
            double T_triple_star = 0.0;
        } troe;

        bool isDuplicate = false;
        std::map<std::string, double> orders;

        // PLOG支持
        struct PLOGPoint {
            double pressure;  // atm
            double A;
            double b;
            double Ea;
        };

        bool isPressureDependent = false;
        std::vector<PLOGPoint> plogData;
    };

    // 热力学数据结构
    struct ThermoData {
        std::string name;
        std::map<std::string, double> composition;

        std::string model;
        std::vector<double> temperatureRanges;

        // NASA7多项式系数
        struct {
            std::vector<double> low;
            std::vector<double> high;
        } coefficients;

        // NASA9多项式数据
        struct NASA9Range {
            std::vector<double> temperatureRange;
            std::vector<double> coefficients;
        };
        std::vector<NASA9Range> nasa9Coeffs;
    };

    // 输运数据结构
    struct TransportData {
        std::string name;
        std::string model;
        std::string geometry;
        double diameter = 0.0;
        double wellDepth = 0.0;
        double dipole = 0.0;
        double polarizability = 0.0;
        double rotationalRelaxation = 0.0;
        std::string note;
    };

    // 整体机理数据结构
    struct MechanismData {
        std::vector<ReactionData> reactions;
        std::vector<ThermoData> thermoSpecies;
        std::vector<TransportData> transportSpecies;
    };

    
    static std::vector<ReactionData> extractKinetics(const std::string& yamlFile, bool verbose = false);
    static std::vector<ThermoData> extractThermo(const std::string& yamlFile, bool verbose = false);
    static std::vector<TransportData> extractTransport(const std::string& yamlFile, bool verbose = false);
    static MechanismData loadMechanism(const std::string& yamlFile, bool verbose = false);

    // 分析和打印函数
    static void analyzeKinetics(const std::string& yamlFile);
    static void analyzeThermo(const std::string& yamlFile);
    static void analyzeTransport(const std::string& yamlFile);
    static void analyzeMechanism(const MechanismData& mechanism, bool printDetails = true, int maxReactions = -1);
    static void printMechanismSummary(const std::string& yamlFile, bool printDetails = true, int maxReactions = -1);

    // 工具函数
    static void parseReactionEquation(const std::string& equation,
        std::map<std::string, double>& reactants,
        std::map<std::string, double>& products);

//private:
//    
//    MechanismData m_mechanismData;
};
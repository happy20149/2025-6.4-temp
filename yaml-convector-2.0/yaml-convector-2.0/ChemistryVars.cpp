#include "ChemistryVars.h"
#include "ChemistryIO.h"
#include <iostream>
#include <sstream>

// 解析动力学数据并返回结构化结果
std::vector<ChemistryVars::ReactionData> ChemistryVars::extractKinetics(const std::string& yamlFile, bool verbose) {
    std::vector<ReactionData> results;

    try {
        // 加载YAML文件
        if (verbose) std::cout << "加载化学动力学文件: " << yamlFile << std::endl;
        ChemistryIO::YamlValue doc = ChemistryIO::loadFile(yamlFile);

        if (!doc.isMap()) {
            std::cerr << "错误: YAML根节点必须是映射表类型" << std::endl;
            return results;
        }
        /*在 YAML 文件格式中，根节点通常是一个映射表（map）结构，即键值对集合。
        ChemistryIO::YamlValue 类封装了 YAML 节点的不同类型，而 asMap() 方法将节点转换为 C++ 中
        的 std::map<std::string, YamlValue> 类型，便于后续访问。*/
        const auto& root = doc.asMap();

        // 检查是否存在反应节点
        if (!root.count("reactions")) {
            if (verbose) std::cout << "未找到反应数据" << std::endl;
            return results;
        }

        // 获取反应列表
        const auto& reactions = root.at("reactions").asSequence();
        if (verbose) std::cout << "找到 " << reactions.size() << " 个反应" << std::endl;

        // 遍历所有反应
        for (size_t i = 0; i < reactions.size(); i++) {
            try {
                const auto& reaction = reactions[i];
                if (!reaction.isMap()) continue;

                const auto& rxnData = reaction.asMap();
                ReactionData reactionItem;

                // 反应方程式
                if (rxnData.count("equation")) {
                    try {
                        reactionItem.equation = rxnData.at("equation").asString();
                        if (verbose) std::cout << "  方程式: " << reactionItem.equation << std::endl;
                    }
                    catch (const std::exception& e) {
                        if (verbose) std::cerr << "  方程式错误: " << e.what() << std::endl;

                        //// 处理特殊情况
                        //if (rxnData.at("equation").isNumber()) {
                        //    double numPrefix = rxnData.at("equation").asNumber();
                        //    if (verbose) std::cout << "  (实际是数值类型: " << numPrefix << ")" << std::endl;

                        //    // 尝试重建反应方程式
                        //    int reactionNum = static_cast<int>(i + 1);
                        //    if (reactionNum == 4) {
                        //        reactionItem.equation = "2 O + M <=> O2 + M";
                        //    }
                        //    else if (reactionNum == 134) {
                        //        reactionItem.equation = "2 CH3 <=> H + C2H5";
                        //    }
                        //    else {
                        //        reactionItem.equation = std::to_string(static_cast<int>(numPrefix)) + " [未知反应]";
                        //    }

                        //    if (verbose) std::cout << "  重建方程式: " << reactionItem.equation << std::endl;
                        //}
                    }
                }

                // 反应类型
                if (rxnData.count("type")) {
                    try {
                        reactionItem.type = rxnData.at("type").asString();
                        if (verbose) std::cout << "  类型: " << reactionItem.type << std::endl;
                    }
                    catch (const std::exception&) {
                        if (verbose) std::cerr << "  类型字段格式错误" << std::endl;
                    }
                }

                // 高压极限速率常数 (用于falloff反应)
                if (rxnData.count("high-P-rate-constant") && rxnData.at("high-P-rate-constant").isMap()) {
                    const auto& highPRate = rxnData.at("high-P-rate-constant").asMap();

                    if (verbose) std::cout << "  高压极限速率常数:" << std::endl;

                    if (highPRate.count("A")) {
                        try {
                            reactionItem.rateConstant.A = highPRate.at("A").asNumber();
                            if (verbose) std::cout << "    A = " << reactionItem.rateConstant.A;

                            if (highPRate.count("A-units")) {
                                reactionItem.rateConstant.A_units = highPRate.at("A-units").asString();
                                if (verbose) std::cout << " " << reactionItem.rateConstant.A_units;
                            }

                            if (verbose) std::cout << std::endl;
                        }
                        catch (const std::exception&) {
                            if (verbose) std::cerr << "    A参数格式错误" << std::endl;
                        }
                    }

                    if (highPRate.count("b")) {
                        try {
                            reactionItem.rateConstant.b = highPRate.at("b").asNumber();
                            if (verbose) std::cout << "    b = " << reactionItem.rateConstant.b << std::endl;
                        }
                        catch (const std::exception&) {
                            if (verbose) std::cerr << "    b参数格式错误" << std::endl;
                        }
                    }

                    if (highPRate.count("Ea")) {
                        try {
                            reactionItem.rateConstant.Ea = highPRate.at("Ea").asNumber();
                            if (verbose) std::cout << "    Ea = " << reactionItem.rateConstant.Ea;

                            if (highPRate.count("Ea-units")) {
                                reactionItem.rateConstant.Ea_units = highPRate.at("Ea-units").asString();
                                if (verbose) std::cout << " " << reactionItem.rateConstant.Ea_units;
                            }

                            if (verbose) std::cout << std::endl;
                        }
                        catch (const std::exception&) {
                            if (verbose) std::cerr << "    Ea参数格式错误" << std::endl;
                        }
                    }
                }


                // 普通反应的速率常数或falloff反应的高压极限（如果high-P-rate-constant未提供）
                else if (rxnData.count("rate-constant") && rxnData.at("rate-constant").isMap()) {
                    const auto& rate = rxnData.at("rate-constant").asMap();

                    if (verbose) {
                        if (reactionItem.type == "falloff" || reactionItem.type == "chemically-activated" || reactionItem.type == "three-body") {
                            std::cout << "  高压极限速率常数 (从rate-constant读取):" << std::endl;
                        }
                        else {
                            std::cout << "  速率常数:" << std::endl;
                        }
                    }

                    if (rate.count("A")) {
                        try {
                            reactionItem.rateConstant.A = rate.at("A").asNumber();
                            if (verbose) std::cout << "    A = " << reactionItem.rateConstant.A;

                            if (rate.count("A-units")) {
                                reactionItem.rateConstant.A_units = rate.at("A-units").asString();
                                if (verbose) std::cout << " " << reactionItem.rateConstant.A_units;
                            }

                            if (verbose) std::cout << std::endl;
                        }
                        catch (const std::exception&) {
                            if (verbose) std::cerr << "    A参数格式错误" << std::endl;
                        }
                    }

                    if (rate.count("b")) {
                        try {
                            reactionItem.rateConstant.b = rate.at("b").asNumber();
                            if (verbose) std::cout << "    b = " << reactionItem.rateConstant.b << std::endl;
                        }
                        catch (const std::exception&) {
                            if (verbose) std::cerr << "    b参数格式错误" << std::endl;
                        }
                    }

                    if (rate.count("Ea")) {
                        try {
                            reactionItem.rateConstant.Ea = rate.at("Ea").asNumber();
                            if (verbose) std::cout << "    Ea = " << reactionItem.rateConstant.Ea;

                            if (rate.count("Ea-units")) {
                                reactionItem.rateConstant.Ea_units = rate.at("Ea-units").asString();
                                if (verbose) std::cout << " " << reactionItem.rateConstant.Ea_units;
                            }

                            if (verbose) std::cout << std::endl;
                        }
                        catch (const std::exception&) {
                            if (verbose) std::cerr << "    Ea参数格式错误" << std::endl;
                        }
                    }
                }

                // 第三体效应
                if (rxnData.count("efficiencies") && rxnData.at("efficiencies").isMap()) {
                    const auto& effs = rxnData.at("efficiencies").asMap();
                    if (verbose) std::cout << "  第三体效率:" << std::endl;

                    //for (const auto& [species, eff] : effs) {
                    for (const auto& pair : effs) {
                        const auto& species = pair.first;
                        const auto& eff = pair.second;
                        try {
                            double value = eff.asNumber();
                            reactionItem.efficiencies[species] = value;
                            if (verbose) std::cout << "    " << species << ": " << value << std::endl;
                        }
                        catch (const std::exception&) {
                            if (verbose) std::cerr << "    " << species << ": 格式错误" << std::endl;
                        }
                    }
                }

                // 低压极限
                if (rxnData.count("low-P-rate-constant") && rxnData.at("low-P-rate-constant").isMap()) {
                    const auto& lowP = rxnData.at("low-P-rate-constant").asMap();
                    if (verbose) std::cout << "  低压极限速率常数:" << std::endl;

                    if (lowP.count("A")) {
                        try {
                            reactionItem.lowPressure.A = lowP.at("A").asNumber();
                            if (verbose) std::cout << "    A = " << reactionItem.lowPressure.A << std::endl;
                        }
                        catch (const std::exception&) {
                            if (verbose) std::cerr << "    A参数格式错误" << std::endl;
                        }
                    }

                    if (lowP.count("b")) {
                        try {
                            reactionItem.lowPressure.b = lowP.at("b").asNumber();
                            if (verbose) std::cout << "    b = " << reactionItem.lowPressure.b << std::endl;
                        }
                        catch (const std::exception&) {
                            if (verbose) std::cerr << "    b参数格式错误" << std::endl;
                        }
                    }

                    if (lowP.count("Ea")) {
                        try {
                            reactionItem.lowPressure.Ea = lowP.at("Ea").asNumber();
                            if (verbose) std::cout << "    Ea = " << reactionItem.lowPressure.Ea << std::endl;
                        }
                        catch (const std::exception&) {
                            if (verbose) std::cerr << "    Ea参数格式错误" << std::endl;
                        }
                    }
                }

                // Troe参数 - 使用正确的字段名
                if (rxnData.count("Troe") && rxnData.at("Troe").isMap()) {
                    const auto& troe = rxnData.at("Troe").asMap();
                    if (verbose) std::cout << "  Troe参数:" << std::endl;

                    // 处理 A 参数
                    if (troe.count("A")) {
                        try {
                            reactionItem.troe.a = troe.at("A").asNumber();
                            if (verbose) std::cout << "    a = " << reactionItem.troe.a << std::endl;
                        }
                        catch (const std::exception&) {
                            if (verbose) std::cerr << "    A参数格式错误" << std::endl;
                        }
                    }

                    // 处理 T3 参数 (对应 T*)
                    if (troe.count("T3")) {
                        try {
                            reactionItem.troe.T_star = troe.at("T3").asNumber();
                            if (verbose) std::cout << "    T* = " << reactionItem.troe.T_star << std::endl;
                        }
                        catch (const std::exception&) {
                            if (verbose) std::cerr << "    T3参数格式错误" << std::endl;
                        }
                    }

                    // 处理 T1 参数 (对应 T**)
                    if (troe.count("T1")) {
                        try {
                            reactionItem.troe.T_double_star = troe.at("T1").asNumber();
                            if (verbose) std::cout << "    T** = " << reactionItem.troe.T_double_star << std::endl;
                        }
                        catch (const std::exception&) {
                            if (verbose) std::cerr << "    T1参数格式错误" << std::endl;
                        }
                    }

                    // 处理 T2 参数 (对应 T***)
                    if (troe.count("T2")) {
                        try {
                            reactionItem.troe.T_triple_star = troe.at("T2").asNumber();
                            if (verbose) std::cout << "    T*** = " << reactionItem.troe.T_triple_star << std::endl;
                        }
                        catch (const std::exception&) {
                            if (verbose) std::cerr << "    T2参数格式错误" << std::endl;
                        }
                    }

                    // 兼容性检查 - 同时支持原始字段名
                    // T*** 参数
                    if (troe.count("T***")) {
                        try {
                            reactionItem.troe.T_triple_star = troe.at("T***").asNumber();
                            if (verbose) std::cout << "    T*** = " << reactionItem.troe.T_triple_star << " (从T***字段读取)" << std::endl;
                        }
                        catch (const std::exception&) {
                            if (verbose) std::cerr << "    T***参数格式错误" << std::endl;
                        }
                    }

                    // T* 参数
                    if (troe.count("T*")) {
                        try {
                            reactionItem.troe.T_star = troe.at("T*").asNumber();
                            if (verbose) std::cout << "    T* = " << reactionItem.troe.T_star << " (从T*字段读取)" << std::endl;
                        }
                        catch (const std::exception&) {
                            if (verbose) std::cerr << "    T*参数格式错误" << std::endl;
                        }
                    }

                    // T** 参数
                    if (troe.count("T**")) {
                        try {
                            reactionItem.troe.T_double_star = troe.at("T**").asNumber();
                            if (verbose) std::cout << "    T** = " << reactionItem.troe.T_double_star << " (从T**字段读取)" << std::endl;
                        }
                        catch (const std::exception&) {
                            if (verbose) std::cerr << "    T**参数格式错误" << std::endl;
                        }
                    }

                    // 兼容性检查 - a字段 vs A字段
                    if (troe.count("a")) {
                        try {
                            reactionItem.troe.a = troe.at("a").asNumber();
                            if (verbose) std::cout << "    a = " << reactionItem.troe.a << " (从a字段读取)" << std::endl;
                        }
                        catch (const std::exception&) {
                            if (verbose) std::cerr << "    a参数格式错误" << std::endl;
                        }
                    }
                }

                // 复制反应
                reactionItem.isDuplicate = rxnData.count("duplicate");
                if (reactionItem.isDuplicate && verbose) {
                    std::cout << "  复制反应: 是" << std::endl;
                }

                // 特殊反应级数
                if (rxnData.count("orders") && rxnData.at("orders").isMap()) {
                    const auto& orders = rxnData.at("orders").asMap();
                    if (verbose) std::cout << "  特殊反应级数:" << std::endl;

                    //for (const auto& [species, order] : orders) {
                    for (const auto& pair : orders) {
                        const auto& species = pair.first;
                        const auto& order = pair.second;
                        try {
                            double value = order.asNumber();
                            reactionItem.orders[species] = value;
                            if (verbose) std::cout << "    " << species << ": " << value << std::endl;
                        }
                        catch (const std::exception&) {
                            if (verbose) std::cerr << "    " << species << ": 格式错误" << std::endl;
                        }
                    }
                }

                // 处理PLOG数据
                if (rxnData.count("note") && rxnData.at("note").isString()) {
                    std::string note = rxnData.at("note").asString();

                    // 检查note中是否包含PLOG数据
                    if (note.find("PLOG/") != std::string::npos) {
                        if (verbose) std::cout << "  发现PLOG压力依赖反应参数" << std::endl;

                        // 将note按行分割
                        std::stringstream ss(note);
                        std::string line;
                        while (std::getline(ss, line)) {
                            // 如果行包含PLOG
                            if (line.find("PLOG/") != std::string::npos) {
                                // 解析PLOG行
                                std::string dataStr = line.substr(line.find("PLOG/") + 5);
                                // 去除末尾的斜杠
                                if (dataStr.back() == '/') {
                                    dataStr.pop_back();
                                }

                                // 分割参数
                                std::stringstream dataStream(dataStr);
                                double pressure, A, b, Ea;
                                if (dataStream >> pressure >> A >> b >> Ea) {
                                    if (verbose) std::cout << "    压力: " << pressure << " atm, A: " << A
                                        << ", b: " << b << ", Ea: " << Ea << std::endl;

                                    // 添加到PLOG数据
                                    reactionItem.rateConstant.isPressureDependent = true;
                                    reactionItem.rateConstant.plogData.push_back({ pressure, A, b, Ea });
                                }
                            }
                        }

                        // 如果找到PLOG数据，则将反应类型更新为pressure-dependent-Arrhenius
                        if (reactionItem.rateConstant.isPressureDependent) {
                            reactionItem.type = "pressure-dependent-Arrhenius";
                            if (verbose) std::cout << "  反应类型更新为: " << reactionItem.type << std::endl;
                        }
                    }
                }

                // 添加到结果集
                results.push_back(reactionItem);
            }
            catch (const std::exception& e) {
                if (verbose) {
                    std::cerr << "处理反应 #" << (i + 1) << " 时出错: " << e.what() << std::endl;
                    std::cerr << "继续处理下一个反应..." << std::endl;
                }
            }
        }
    }
    catch (const std::exception& e) {
        std::cerr << "错误: " << e.what() << std::endl;
    }

    return results;
}

// 解析热力学数据并返回结构化结果
std::vector<ChemistryVars::ThermoData> ChemistryVars::extractThermo(const std::string& yamlFile, bool verbose) {
    std::vector<ThermoData> results;

    try {
        // 加载YAML文件
        if (verbose) std::cout << "加载热力学数据文件: " << yamlFile << std::endl;
        ChemistryIO::YamlValue doc = ChemistryIO::loadFile(yamlFile);

        if (!doc.isMap()) {
            std::cerr << "错误: YAML根节点必须是映射表类型" << std::endl;
            return results;
        }

        const auto& root = doc.asMap();

        // 检查是否存在组分节点
        if (!root.count("species")) {
            if (verbose) std::cout << "未找到组分数据" << std::endl;
            return results;
        }

        // 获取组分列表
        const auto& speciesList = root.at("species").asSequence();
        if (verbose) std::cout << "找到 " << speciesList.size() << " 个组分" << std::endl;

        // 统计有效NASA7数据的组分数量
        int validNASA7Count = 0;

        // 遍历所有组分
        for (size_t i = 0; i < speciesList.size(); i++) {
            try {
                const auto& species = speciesList[i];
                if (!species.isMap()) continue;

                const auto& speciesData = species.asMap();
                ThermoData thermoItem;

                if (verbose) std::cout << "\n组分 #" << (i + 1) << ":" << std::endl;

                // 组分名称
                if (speciesData.count("name")) {
                    try {
                        thermoItem.name = speciesData.at("name").asString();
                        if (verbose) std::cout << "  名称: " << thermoItem.name << std::endl;
                    }
                    catch (const std::exception&) {
                        if (verbose) std::cerr << "  名称格式错误" << std::endl;
                    }
                }

                // 组分组成
                if (speciesData.count("composition") && speciesData.at("composition").isMap()) {
                    const auto& composition = speciesData.at("composition").asMap();
                    if (verbose) std::cout << "  组成: ";

                    //for (const auto& [element, count] : composition) {
                    for (const auto& pair : composition) {
                        const auto& element = pair.first;
                        const auto& count = pair.second;
                        try {
                            double value = count.asNumber();
                            thermoItem.composition[element] = value;
                            if (verbose) std::cout << element << ":" << value << " ";
                        }
                        catch (const std::exception&) {
                            if (verbose) std::cout << element << ":[格式错误] ";
                        }
                    }

                    if (verbose) std::cout << std::endl;
                }

                // 热力学数据
                bool hasNASA7Data = false;
                if (speciesData.count("thermo") && speciesData.at("thermo").isMap()) {
                    const auto& thermo = speciesData.at("thermo").asMap();
                    if (verbose) std::cout << "  热力学数据:" << std::endl;

                    // 热力学模型
                    if (thermo.count("model")) {
                        try {
                            thermoItem.model = thermo.at("model").asString();
                            if (verbose) std::cout << "    模型: " << thermoItem.model << std::endl;

                            // 检查是否为NASA7模型
                            hasNASA7Data = (thermoItem.model == "NASA7");
                        }
                        catch (const std::exception&) {
                            if (verbose) std::cerr << "    模型格式错误" << std::endl;
                        }
                    }

                    // 温度范围
                    bool validTempRanges = false;
                    if (thermo.count("temperature-ranges") && thermo.at("temperature-ranges").isSequence()) {
                        const auto& tempRanges = thermo.at("temperature-ranges").asSequence();

                        // 清空已有数据
                        thermoItem.temperatureRanges.clear();

                        if (verbose) std::cout << "    温度范围(K): ";

                        for (const auto& temp : tempRanges) {
                            try {
                                double value = temp.asNumber();
                                thermoItem.temperatureRanges.push_back(value);
                                if (verbose) std::cout << value << " ";
                            }
                            catch (const std::exception&) {
                                if (verbose) std::cout << "[格式错误] ";
                            }
                        }

                        if (verbose) std::cout << std::endl;

                        // 验证NASA7温度范围（应该有3个点，对应2个温度区间）
                        if (hasNASA7Data && thermoItem.temperatureRanges.size() == 3) {
                            validTempRanges = true;
                            if (verbose) std::cout << "    NASA7温度范围有效: "
                                << thermoItem.temperatureRanges[0] << " - "
                                << thermoItem.temperatureRanges[1] << " - "
                                << thermoItem.temperatureRanges[2] << " K" << std::endl;
                        }
                        else if (hasNASA7Data) {
                            if (verbose) std::cerr << "    警告: NASA7模型需要3个温度点，实际有 "
                                << thermoItem.temperatureRanges.size() << std::endl;
                        }
                    }
                    else if (hasNASA7Data && verbose) {
                        std::cerr << "    警告: 缺少NASA7温度范围数据" << std::endl;
                    }

                    // NASA多项式系数
                    bool validCoeffs = false;
                    bool validLowCoeffs = false;
                    bool validHighCoeffs = false;

                    // 检查标准格式 (coefficients.low/high)
                    if (thermo.count("coefficients") && thermo.at("coefficients").isMap()) {
                        const auto& coeffs = thermo.at("coefficients").asMap();
                        if (verbose) std::cout << "    NASA多项式系数:" << std::endl;

                        // 清空已有系数数据
                        thermoItem.coefficients.low.clear();
                        thermoItem.coefficients.high.clear();

                        // 低温系数
                        if (coeffs.count("low") && coeffs.at("low").isSequence()) {
                            const auto& lowCoeffs = coeffs.at("low").asSequence();
                            if (verbose) std::cout << "      低温系数: ";

                            for (const auto& coeff : lowCoeffs) {
                                try {
                                    double value = coeff.asNumber();
                                    thermoItem.coefficients.low.push_back(value);
                                    if (verbose) std::cout << value << " ";
                                }
                                catch (const std::exception&) {
                                    if (verbose) std::cout << "[格式错误] ";
                                }
                            }

                            if (verbose) std::cout << std::endl;

                            // 验证NASA7低温系数（应该是7个）
                            if (hasNASA7Data && thermoItem.coefficients.low.size() == 7) {
                                validLowCoeffs = true;
                                if (verbose) std::cout << "      NASA7低温系数有效" << std::endl;
                            }
                            else if (hasNASA7Data) {
                                if (verbose) std::cerr << "      警告: NASA7模型需要7个低温系数，实际有 "
                                    << thermoItem.coefficients.low.size() << std::endl;
                            }
                        }
                        else if (hasNASA7Data && verbose) {
                            std::cerr << "      警告: 缺少NASA7低温系数数据" << std::endl;
                        }

                        // 高温系数
                        if (coeffs.count("high") && coeffs.at("high").isSequence()) {
                            const auto& highCoeffs = coeffs.at("high").asSequence();
                            if (verbose) std::cout << "      高温系数: ";

                            for (const auto& coeff : highCoeffs) {
                                try {
                                    double value = coeff.asNumber();
                                    thermoItem.coefficients.high.push_back(value);
                                    if (verbose) std::cout << value << " ";
                                }
                                catch (const std::exception&) {
                                    if (verbose) std::cout << "[格式错误] ";
                                }
                            }

                            if (verbose) std::cout << std::endl;

                            // 验证NASA7高温系数（应该是7个）
                            if (hasNASA7Data && thermoItem.coefficients.high.size() == 7) {
                                validHighCoeffs = true;
                                if (verbose) std::cout << "      NASA7高温系数有效" << std::endl;
                            }
                            else if (hasNASA7Data) {
                                if (verbose) std::cerr << "      警告: NASA7模型需要7个高温系数，实际有 "
                                    << thermoItem.coefficients.high.size() << std::endl;
                            }
                        }
                        else if (hasNASA7Data && verbose) {
                            std::cerr << "      警告: 缺少NASA7高温系数数据" << std::endl;
                        }
                    }
                    // 检查替代格式 (data 数组格式)
                    else if (thermo.count("data") && thermo.at("data").isSequence()) {
                        const auto& dataArray = thermo.at("data").asSequence();
                        if (verbose) std::cout << "    NASA多项式系数 (data格式):" << std::endl;

                        // 清空已有系数数据
                        thermoItem.coefficients.low.clear();
                        thermoItem.coefficients.high.clear();

                        // 读取低温系数 (第一组数据)
                        if (dataArray.size() > 0 && dataArray[0].isSequence()) {
                            const auto& lowCoeffs = dataArray[0].asSequence();
                            if (verbose) std::cout << "      低温系数: ";

                            for (const auto& coeff : lowCoeffs) {
                                try {
                                    double value = coeff.asNumber();
                                    thermoItem.coefficients.low.push_back(value);
                                    if (verbose) std::cout << value << " ";
                                }
                                catch (const std::exception&) {
                                    if (verbose) std::cout << "[格式错误] ";
                                }
                            }
                            if (verbose) std::cout << std::endl;

                            // 验证NASA7低温系数
                            if (hasNASA7Data && thermoItem.coefficients.low.size() == 7) {
                                validLowCoeffs = true;
                                if (verbose) std::cout << "      NASA7低温系数有效" << std::endl;
                            }
                            else if (hasNASA7Data) {
                                if (verbose) std::cerr << "      警告: NASA7模型需要7个低温系数，实际有 "
                                    << thermoItem.coefficients.low.size() << std::endl;
                            }
                        }
                        else if (hasNASA7Data && verbose) {
                            std::cerr << "      警告: 缺少NASA7低温系数数据" << std::endl;
                        }

                        // 读取高温系数 (第二组数据)
                        if (dataArray.size() > 1 && dataArray[1].isSequence()) {
                            const auto& highCoeffs = dataArray[1].asSequence();
                            if (verbose) std::cout << "      高温系数: ";

                            for (const auto& coeff : highCoeffs) {
                                try {
                                    double value = coeff.asNumber();
                                    thermoItem.coefficients.high.push_back(value);
                                    if (verbose) std::cout << value << " ";
                                }
                                catch (const std::exception&) {
                                    if (verbose) std::cout << "[格式错误] ";
                                }
                            }
                            if (verbose) std::cout << std::endl;

                            // 验证NASA7高温系数
                            if (hasNASA7Data && thermoItem.coefficients.high.size() == 7) {
                                validHighCoeffs = true;
                                if (verbose) std::cout << "      NASA7高温系数有效" << std::endl;
                            }
                            else if (hasNASA7Data) {
                                if (verbose) std::cerr << "      警告: NASA7模型需要7个高温系数，实际有 "
                                    << thermoItem.coefficients.high.size() << std::endl;
                            }
                        }
                        else if (hasNASA7Data && verbose) {
                            std::cerr << "      警告: 缺少NASA7高温系数数据" << std::endl;
                        }
                    }
                    else if (hasNASA7Data && verbose) {
                        std::cerr << "    警告: 缺少NASA7系数数据 (未找到 'coefficients' 或 'data' 字段)" << std::endl;
                    }

                    // 检查系数是否完整有效
                    validCoeffs = validLowCoeffs && validHighCoeffs;

                    // 汇总NASA7数据有效性
                    if (hasNASA7Data && validTempRanges && validCoeffs) {
                        validNASA7Count++;
                        if (verbose) std::cout << "    NASA7数据完整有效" << std::endl;
                    }
                }

                // NASA-9多项式格式支持
                if (speciesData.count("nasa9-coeffs") && speciesData.at("nasa9-coeffs").isSequence()) {
                    const auto& nasa9Ranges = speciesData.at("nasa9-coeffs").asSequence();
                    if (verbose) std::cout << "  NASA-9多项式数据:" << std::endl;

                    // 清空NASA9系数数据
                    thermoItem.nasa9Coeffs.clear();

                    for (size_t j = 0; j < nasa9Ranges.size(); j++) {
                        try {
                            const auto& range = nasa9Ranges[j].asMap();
                            ThermoData::NASA9Range nasa9Range;

                            if (verbose) std::cout << "    温度范围 #" << (j + 1) << ":" << std::endl;

                            if (range.count("T-range")) {
                                try {
                                    const auto& tRange = range.at("T-range").asSequence();
                                    double tMin = tRange[0].asNumber();
                                    double tMax = tRange[1].asNumber();

                                    nasa9Range.temperatureRange.push_back(tMin);
                                    nasa9Range.temperatureRange.push_back(tMax);

                                    if (verbose) std::cout << "      温度: " << tMin << " - " << tMax << " K" << std::endl;
                                }
                                catch (const std::exception&) {
                                    if (verbose) std::cerr << "      温度范围格式错误" << std::endl;
                                }
                            }

                            if (range.count("coeffs")) {
                                try {
                                    const auto& rangeCoeffs = range.at("coeffs").asSequence();
                                    if (verbose) std::cout << "      系数: ";

                                    for (const auto& coeff : rangeCoeffs) {
                                        double value = coeff.asNumber();
                                        nasa9Range.coefficients.push_back(value);
                                        if (verbose) std::cout << value << " ";
                                    }

                                    if (verbose) std::cout << std::endl;

                                    // 验证NASA9系数数量（应为9个）
                                    if (nasa9Range.coefficients.size() != 9 && verbose) {
                                        std::cerr << "      警告: NASA9模型需要9个系数，实际有 "
                                            << nasa9Range.coefficients.size() << std::endl;
                                    }
                                }
                                catch (const std::exception&) {
                                    if (verbose) std::cerr << "      系数格式错误" << std::endl;
                                }
                            }

                            thermoItem.nasa9Coeffs.push_back(nasa9Range);
                        }
                        catch (const std::exception&) {
                            if (verbose) std::cerr << "    处理NASA9温度范围 #" << (j + 1) << " 时出错" << std::endl;
                        }
                    }
                }

                // 添加到结果集
                results.push_back(thermoItem);
            }
            catch (const std::exception& e) {
                if (verbose) {
                    std::cerr << "处理组分 #" << (i + 1) << " 时出错: " << e.what() << std::endl;
                    std::cerr << "继续处理下一个组分..." << std::endl;
                }
            }
        }

        // 输出NASA7数据统计
        if (verbose) {
            std::cout << "\n总结: 找到 " << results.size() << " 个组分，其中 "
                << validNASA7Count << " 个有完整有效的NASA7热力学数据" << std::endl;
        }
    }
    catch (const std::exception& e) {
        std::cerr << "错误: " << e.what() << std::endl;
    }

    return results;
}

// 解析输运性质数据并返回结构化结果
std::vector<ChemistryVars::TransportData> ChemistryVars::extractTransport(const std::string& yamlFile, bool verbose) {
    std::vector<TransportData> results;

    try {
        // 加载YAML文件
        if (verbose) std::cout << "加载输运性质数据文件: " << yamlFile << std::endl;
        ChemistryIO::YamlValue doc = ChemistryIO::loadFile(yamlFile);

        if (!doc.isMap()) {
            std::cerr << "错误: YAML根节点必须是映射表类型" << std::endl;
            return results;
        }

        const auto& root = doc.asMap();

        // 检查是否存在组分节点
        if (!root.count("species")) {
            if (verbose) std::cout << "未找到组分数据" << std::endl;
            return results;
        }

        // 获取组分列表
        const auto& speciesList = root.at("species").asSequence();
        if (verbose) std::cout << "找到 " << speciesList.size() << " 个组分" << std::endl;

        int speciesWithTransport = 0;

        // 遍历所有组分
        for (size_t i = 0; i < speciesList.size(); i++) {
            try {
                const auto& species = speciesList[i];
                if (!species.isMap()) continue;

                const auto& speciesData = species.asMap();

                // 仅处理有输运数据的组分
                if (!speciesData.count("transport") || !speciesData.at("transport").isMap()) {
                    continue;
                }

                speciesWithTransport++;
                TransportData transportItem;

                // 组分名称
                if (speciesData.count("name")) {
                    try {
                        transportItem.name = speciesData.at("name").asString();
                    }
                    catch (const std::exception&) {
                        transportItem.name = "未知组分";
                    }
                }
                else {
                    transportItem.name = "未知组分";
                }

                if (verbose) {
                    std::cout << "\n组分 #" << (i + 1) << " (" << transportItem.name << ") 输运性质:" << std::endl;
                }

                // 获取输运数据
                const auto& transport = speciesData.at("transport").asMap();

                // 输运模型
                if (transport.count("model")) {
                    try {
                        transportItem.model = transport.at("model").asString();
                        if (verbose) std::cout << "  模型: " << transportItem.model << std::endl;
                    }
                    catch (const std::exception&) {
                        if (verbose) std::cerr << "  模型格式错误" << std::endl;
                    }
                }

                // 几何构型
                if (transport.count("geometry")) {
                    try {
                        transportItem.geometry = transport.at("geometry").asString();
                        if (verbose) std::cout << "  几何构型: " << transportItem.geometry << std::endl;
                    }
                    catch (const std::exception&) {
                        if (verbose) std::cerr << "  几何构型格式错误" << std::endl;
                    }
                }

                // 碰撞直径
                if (transport.count("diameter")) {
                    try {
                        transportItem.diameter = transport.at("diameter").asNumber();
                        if (verbose) std::cout << "  碰撞直径: " << transportItem.diameter << " Å" << std::endl;
                    }
                    catch (const std::exception&) {
                        if (verbose) std::cerr << "  碰撞直径格式错误" << std::endl;
                    }
                }

                // 势阱深度
                if (transport.count("well-depth")) {
                    try {
                        transportItem.wellDepth = transport.at("well-depth").asNumber();
                        if (verbose) std::cout << "  势阱深度: " << transportItem.wellDepth << " K" << std::endl;
                    }
                    catch (const std::exception&) {
                        if (verbose) std::cerr << "  势阱深度格式错误" << std::endl;
                    }
                }

                // 偶极矩
                if (transport.count("dipole")) {
                    try {
                        transportItem.dipole = transport.at("dipole").asNumber();
                        if (verbose) std::cout << "  偶极矩: " << transportItem.dipole << " Debye" << std::endl;
                    }
                    catch (const std::exception&) {
                        if (verbose) std::cerr << "  偶极矩格式错误" << std::endl;
                    }
                }

                // 极化率
                if (transport.count("polarizability")) {
                    try {
                        transportItem.polarizability = transport.at("polarizability").asNumber();
                        if (verbose) std::cout << "  极化率: " << transportItem.polarizability << " Å³" << std::endl;
                    }
                    catch (const std::exception&) {
                        if (verbose) std::cerr << "  极化率格式错误" << std::endl;
                    }
                }

                // 转动松弛数
                if (transport.count("rotational-relaxation")) {
                    try {
                        transportItem.rotationalRelaxation = transport.at("rotational-relaxation").asNumber();
                        if (verbose) std::cout << "  转动松弛数: " << transportItem.rotationalRelaxation << std::endl;
                    }
                    catch (const std::exception&) {
                        if (verbose) std::cerr << "  转动松弛数格式错误" << std::endl;
                    }
                }

                // 附加说明
                if (transport.count("note")) {
                    try {
                        transportItem.note = transport.at("note").asString();
                        if (verbose) std::cout << "  附加说明: " << transportItem.note << std::endl;
                    }
                    catch (const std::exception&) {
                        if (verbose) std::cerr << "  附加说明格式错误" << std::endl;
                    }
                }

                // 添加到结果集
                results.push_back(transportItem);

            }
            catch (const std::exception& e) {
                if (verbose) {
                    std::cerr << "处理组分 #" << (i + 1) << " 输运性质时出错: " << e.what() << std::endl;
                    std::cerr << "继续处理下一个组分..." << std::endl;
                }
            }
        }

        if (verbose) {
            std::cout << "\n总计: " << speciesWithTransport << " 个组分具有输运性质数据" << std::endl;
        }

    }
    catch (const std::exception& e) {
        std::cerr << "错误: " << e.what() << std::endl;
    }

    return results;
}

// 加载整个机理数据
ChemistryVars::MechanismData ChemistryVars::loadMechanism(const std::string& yamlFile, bool verbose) {
    MechanismData mechanism;

    mechanism.reactions = extractKinetics(yamlFile, verbose);
    mechanism.thermoSpecies = extractThermo(yamlFile, verbose);
    mechanism.transportSpecies = extractTransport(yamlFile, verbose);

    return mechanism;
}

// 保留原有的分析函数 - 直接调用extract函数并显示
void ChemistryVars::analyzeKinetics(const std::string& yamlFile) {
    extractKinetics(yamlFile, true);
}

void ChemistryVars::analyzeThermo(const std::string& yamlFile) {
    extractThermo(yamlFile, true);
}

void ChemistryVars::analyzeTransport(const std::string& yamlFile) {
    extractTransport(yamlFile, true);
}

// 解析反应方程式，提取反应物、产物及其化学计量数
void ChemistryVars::parseReactionEquation(const std::string& equation,
    std::map<std::string, double>& reactants,
    std::map<std::string, double>& products) {
    // 清空输入映射
    reactants.clear();
    products.clear();

    // 定位反应箭头
    size_t arrowPos = equation.find("<=>");
    if (arrowPos == std::string::npos) {
        arrowPos = equation.find("=>");
        if (arrowPos == std::string::npos) {
            arrowPos = equation.find("=");
        }
    }

    if (arrowPos == std::string::npos) return; // 未找到反应箭头

    // 分割反应物和产物
    std::string reactantsStr = equation.substr(0, arrowPos);
    std::string productsStr = equation.substr(arrowPos +
        (equation.substr(arrowPos).find(">") != std::string::npos ? 3 : 1));

    // 解析函数
    auto parseSpecies = [](const std::string& side, std::map<std::string, double>& species) {
        std::stringstream ss(side);
        std::string token;
        double stoich = 1.0;
        bool expectSpecies = true;

        while (ss >> token) {
            if (token == "+") {
                expectSpecies = true;
                continue;
            }

            // 检查是否是数字开头
            if (isdigit(token[0])) {
                size_t endPos;
                stoich = std::stod(token, &endPos);

                // 如果整个token是数字，等待下一个token作为组分名
                if (endPos == token.length()) {
                    expectSpecies = false;
                    continue;
                }

                // 否则，剩余部分是组分名
                std::string speciesName = token.substr(endPos);
                species[speciesName] += stoich;
                stoich = 1.0;
                expectSpecies = true;
            }
            else if (!expectSpecies) {
                // 前面读到了系数，这个是组分名
                species[token] += stoich;
                stoich = 1.0;
                expectSpecies = true;
            }
            else {
                // 没有前导系数，默认为1
                species[token] += 1.0;
            }
        }
        };

    parseSpecies(reactantsStr, reactants);
    parseSpecies(productsStr, products);
}

// 打印机理数据摘要及详细信息
void ChemistryVars::printMechanismSummary(const std::string& yamlFile, bool printDetails, int maxReactions) {
    // 加载机理数据但不打印详细信息
    MechanismData mechanism = loadMechanism(yamlFile, false);

    // 调用处理函数
    analyzeMechanism(mechanism, printDetails, maxReactions); // 化学反应机制
    analyzeThermo(yamlFile); // 热力学
    analyzeTransport(yamlFile); // 输运性质
}

// 处理和分析已加载的机理数据
void ChemistryVars::analyzeMechanism(const MechanismData& mechanism, bool printDetails, int maxReactions) {
    // 打印基本摘要信息
    std::cout << "成功加载机理数据:" << std::endl;
    std::cout << "  " << mechanism.reactions.size() << " 个反应" << std::endl;
    std::cout << "  " << mechanism.thermoSpecies.size() << " 个组分热力学数据" << std::endl;
    std::cout << "  " << mechanism.transportSpecies.size() << " 个组分输运性质数据" << std::endl;

    // 如果不需要打印详细信息，直接返回
    if (!printDetails) {
        return;
    }

    // 遍历并访问所有反应的数据
    int count = 0;
    for (const auto& reaction : mechanism.reactions) {
        std::cout << "反应: " << reaction.equation << std::endl;

        // 显示反应类型
        if (!reaction.type.empty()) {
            std::cout << "  反应类型: " << reaction.type << std::endl;
        }

        // 解析反应方程式获取化学计量数
        std::map<std::string, double> reactants, products;
        parseReactionEquation(reaction.equation, reactants, products);

        // 打印反应物及其化学计量数
        std::cout << "  反应物:" << std::endl;
        //for (const auto& [species, coeff] : reactants) {
        for (const auto& pair : reactants) {
            const auto& species = pair.first;
            const auto& coeff = pair.second;
            std::cout << "    " << species << ": " << coeff << std::endl;
        }

        // 打印产物及其化学计量数
        std::cout << "  产物:" << std::endl;
        //for (const auto& [species, coeff] : products) {
        for (const auto& pair : products) {
            const auto& species = pair.first;
            const auto& coeff = pair.second;
            std::cout << "    " << species << ": " << coeff << std::endl;
        }

        // 根据反应类型显示不同的速率参数
        if (reaction.type == "falloff" || reaction.type == "chemically-activated") {
            // Falloff/化学活化反应 - 显示高压和低压限制参数
            std::cout << "  高压极限速率参数:" << std::endl;
            std::cout << "    A = " << reaction.rateConstant.A;
            if (!reaction.rateConstant.A_units.empty()) {
                std::cout << " " << reaction.rateConstant.A_units;
            }
            std::cout << std::endl;
            std::cout << "    b = " << reaction.rateConstant.b << std::endl;
            std::cout << "    Ea = " << reaction.rateConstant.Ea;
            if (!reaction.rateConstant.Ea_units.empty()) {
                std::cout << " " << reaction.rateConstant.Ea_units;
            }
            std::cout << std::endl;

            // 显示低压极限参数
            std::cout << "  低压极限速率参数:" << std::endl;
            std::cout << "    A = " << reaction.lowPressure.A;
            if (!reaction.rateConstant.A_units.empty()) {
                std::cout << " " << reaction.rateConstant.A_units;
            }
            std::cout << std::endl;
            std::cout << "    b = " << reaction.lowPressure.b << std::endl;
            std::cout << "    Ea = " << reaction.lowPressure.Ea;
            if (!reaction.rateConstant.Ea_units.empty()) {
                std::cout << " " << reaction.rateConstant.Ea_units;
            }
            std::cout << std::endl;

            // 显示Troe参数（如果有）
            if (reaction.troe.a != 0.0 || reaction.troe.T_star != 0.0 ||
                reaction.troe.T_double_star != 0.0 || reaction.troe.T_triple_star != 0.0) {
                std::cout << "  Troe参数:" << std::endl;
                std::cout << "    a = " << reaction.troe.a << std::endl;
                if (reaction.troe.T_star != 0.0) {
                    std::cout << "    T* = " << reaction.troe.T_star << " K" << std::endl;
                }
                if (reaction.troe.T_double_star != 0.0) {
                    std::cout << "    T** = " << reaction.troe.T_double_star << " K" << std::endl;
                }
                if (reaction.troe.T_triple_star != 0.0) {
                    std::cout << "    T*** = " << reaction.troe.T_triple_star << " K" << std::endl;
                }
            }
        }
        else if (reaction.type == "three-body") {
            // 三体反应 - 只显示一组速率参数（在YAML中是rate-constant字段）
            std::cout << "  反应速率参数:" << std::endl;
            std::cout << "    A = " << reaction.rateConstant.A;
            if (!reaction.rateConstant.A_units.empty()) {
                std::cout << " " << reaction.rateConstant.A_units;
            }
            std::cout << std::endl;
            std::cout << "    b = " << reaction.rateConstant.b << std::endl;
            std::cout << "    Ea = " << reaction.rateConstant.Ea;
            if (!reaction.rateConstant.Ea_units.empty()) {
                std::cout << " " << reaction.rateConstant.Ea_units;
            }
            std::cout << std::endl;
        }
        else {
            // 普通反应 - 显示标准速率参数
            std::cout << "  反应速率参数:" << std::endl;
            std::cout << "    A = " << reaction.rateConstant.A;
            if (!reaction.rateConstant.A_units.empty()) {
                std::cout << " " << reaction.rateConstant.A_units;
            }
            std::cout << std::endl;
            std::cout << "    b = " << reaction.rateConstant.b << std::endl;
            std::cout << "    Ea = " << reaction.rateConstant.Ea;
            if (!reaction.rateConstant.Ea_units.empty()) {
                std::cout << " " << reaction.rateConstant.Ea_units;
            }
            std::cout << std::endl;
        }

        // 显示第三体效率（对于任何可能有第三体效率的反应）
        if (!reaction.efficiencies.empty()) {
            std::cout << "  第三体效率:" << std::endl;
            //for (const auto& [species, eff] : reaction.efficiencies) {
            for (const auto& pair : reaction.efficiencies) {
                const auto& species = pair.first;
                const auto& eff = pair.second;
                std::cout << "    " << species << ": " << eff << std::endl;
            }
        }

        // 如果是复制反应，显示标记
        if (reaction.isDuplicate) {
            std::cout << "  [复制反应]" << std::endl;
        }

        // 如果有特殊反应级数，打印出来
        if (!reaction.orders.empty()) {
            std::cout << "  特殊反应级数:" << std::endl;
            //for (const auto& [species, order] : reaction.orders) {
            for (const auto& pair : reaction.orders) {
                const auto& species = pair.first;
                const auto& order = pair.second;
                std::cout << "    " << species << ": " << order << std::endl;
            }
        }

        // 在打印反应信息的部分添加PLOG数据显示
        if (reaction.rateConstant.isPressureDependent) {
            std::cout << "  压力依赖Arrhenius参数(PLOG):" << std::endl;
            for (const auto& plog : reaction.rateConstant.plogData) {
                std::cout << "    压力: " << plog.pressure << " atm" << std::endl;
                std::cout << "      A = " << plog.A;
                if (!reaction.rateConstant.A_units.empty()) {
                    std::cout << " " << reaction.rateConstant.A_units;
                }
                std::cout << std::endl;
                std::cout << "      b = " << plog.b << std::endl;
                std::cout << "      Ea = " << plog.Ea;
                if (!reaction.rateConstant.Ea_units.empty()) {
                    std::cout << " " << reaction.rateConstant.Ea_units;
                }
                std::cout << std::endl;
            }
        }

        std::cout << std::endl;

        // 计数器增加
        count++;

        // 如果达到最大显示数量，跳出循环
        if (maxReactions > 0 && count >= maxReactions) {
            std::cout << "... 仅显示前 " << maxReactions << " 个反应，共 "
                << mechanism.reactions.size() << " 个反应" << std::endl;
            break;
        }
    }
}
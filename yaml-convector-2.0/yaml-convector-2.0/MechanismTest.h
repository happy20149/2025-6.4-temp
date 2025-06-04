#pragma once
#include <string>
#include "ChemistryVars.h"
#include "ChemistryIO.h"

// 声明测试函数（定义在MechanismTest.cpp中）
bool runMechanismTest(const std::string& yamlFile);//测试ch4机理
bool testSpecificC2H4Data(const std::string& yamlFile); // 测试c2h4机理

bool chemkinIOTest(const std::string& yamlFile);//整体测试接口
#include "ChemistryVars.h"
#include "ChemistryIO.h"
#include "IdealGasPhase.h"
#include <iostream>
#include <vector>

int main() {
    try {
        std::cout << "=== Debug Species Loading ===" << std::endl;
        
        // Test ChemistryVars::extractThermo directly
        std::cout << "\n1. Testing ChemistryVars::extractThermo directly:" << std::endl;
        std::vector<ChemistryVars::ThermoData> thermoData = ChemistryVars::extractThermo("D:\\cantera\\data\\h2o2.yaml", true);
        std::cout << "extractThermo returned " << thermoData.size() << " species" << std::endl;
        
        for (size_t i = 0; i < thermoData.size(); ++i) {
            std::cout << "  Species " << (i+1) << ": " << thermoData[i].name << std::endl;
        }
        
        // Test IdealGasPhase initialization
        std::cout << "\n2. Testing IdealGasPhase initialization:" << std::endl;
        YamlConvector2::IdealGasPhase gas;
        gas.initFromYaml("D:\\cantera\\data\\h2o2.yaml", "ohmech");
        
        std::cout << "IdealGasPhase loaded " << gas.nSpecies() << " species" << std::endl;
        for (size_t i = 0; i < gas.nSpecies(); ++i) {
            std::cout << "  Species " << (i+1) << ": " << gas.speciesName(i) << std::endl;
        }
        
        // Test setting state and reporting
        std::cout << "\n3. Testing state setting and reporting:" << std::endl;
        gas.setState_TPX(1500.0, 2.0 * YamlConvector2::OneAtm, "O2:1.0, H2:3.0, AR:1.0");
        
        std::cout << "\n" << gas.report() << std::endl;
        
    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}

#include "ChemistryVars.h"
#include "ChemistryIO.h"
#include "ChemEquil.h"
#include <iostream>
#include <limits>
#include "MechanismTest.h"
#include "IdealGasPhase.h"

void testChemicalEquilibrium() {    std::cout << "=== Testing Chemical Equilibrium ===" << std::endl;
      try {
        YamlConvector2::IdealGasPhase gas;
        gas.initFromYaml("..\\..\\..\\..\\h2o2.yaml", "ohmech");
        
        // Debug: Show how many species were loaded
        std::cout << "Loaded " << gas.nSpecies() << " species from h2o2.yaml:" << std::endl;
        for (size_t i = 0; i < gas.nSpecies(); ++i) {
            std::cout << "  " << (i+1) << ". " << gas.speciesName(i) << std::endl;
        }
        std::cout << std::endl;
        
        // Set initial state - same as cantera_cpp_demo
        gas.setState_TPX(1500.0, 2.0 * YamlConvector2::OneAtm, "O2:1.0, H2:3.0, AR:1.0");
        
        std::cout << "=== Initial State ===" << std::endl;
        std::cout << gas.report() << std::endl;
        
        // Perform equilibrium calculation
        std::cout << "\n=== Performing TP Equilibrium ===" << std::endl;
        int result = gas.equilibrate("TP");
        
        if (result == 0) {
            std::cout << "\n=== Equilibrium State ===" << std::endl;
            std::cout << gas.report() << std::endl;
        } else {
            std::cout << "Equilibrium calculation failed with code: " << result << std::endl;
        }
        
    } catch (const std::exception& e) {
        std::cout << "Error in equilibrium test: " << e.what() << std::endl;
    }
}

int main() {
    testChemicalEquilibrium();
    return 0;
}
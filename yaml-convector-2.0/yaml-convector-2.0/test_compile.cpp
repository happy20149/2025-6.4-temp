#include "ChemEquil.h"
#include "IdealGasPhase.h"
#include <iostream>

int main() {
    std::cout << "Testing ChemEquil compilation..." << std::endl;
    
    try {
        // Test basic instantiation
        YamlConvector2::ChemEquil equilibrium;
        std::cout << "ChemEquil object created successfully" << std::endl;
        
        YamlConvector2::IdealGasPhase gas;
        std::cout << "IdealGasPhase object created successfully" << std::endl;
        
        std::cout << "Basic compilation test passed!" << std::endl;
        return 0;
        
    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
        return 1;
    }
}

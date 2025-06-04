#include "ChemistryVars.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <limits>

// Test results structure
struct TestResults {
    int totalTests = 0;
    int passedTests = 0;
    std::vector<std::string> failureMessages;
};

// Utility function to compare floating-point values
bool isEqual(double a, double b, double tolerance = 1e-5) {
    return std::abs(a - b) < tolerance;
}

// Utility function to compare vectors of doubles
bool compareVectors(const std::vector<double>& actual, const std::vector<double>& expected, double tolerance = 1e-5) {
    if (actual.size() != expected.size()) return false;
    for (size_t i = 0; i < actual.size(); i++) {
        if (!isEqual(actual[i], expected[i], tolerance)) return false;
    }
    return true;
}

// Test thermodynamic data
void testThermo(ChemistryVars::MechanismData& mechanism, TestResults& results) {
    std::cout << "\n===== Thermodynamic Data Tests =====\n";

    // Expected thermodynamic data
    struct ExpectedThermo {
        std::string name;
        std::map<std::string, double> composition;
        std::string model;
        std::vector<double> temperatureRanges;
        std::vector<double> lowCoeffs;
        std::vector<double> highCoeffs;
        std::string note;  // Optional note field
    };

    std::vector<ExpectedThermo> expectedThermo = {
        // H2
        {
            "H2",
            {{"H", 2.0}},
            "NASA7",
            {200.0, 1000.0, 6000.0},
            {2.34433112, 7.98052075e-03, -1.9478151e-05, 2.01572094e-08, -7.37611761e-12, -917.935173, 0.683010238},
            {2.93286575, 8.26608026e-04, -1.46402364e-07, 1.54100414e-11, -6.888048e-16, -813.065581, -1.02432865},
            "TPIS78"
        },
        // // H
        // {
        //     "H",
        //     {{"H", 1.0}},
        //     "NASA7",
        //     {200.0, 6000.0},
        //     {2.5, 0.0, 0.0, 0.0, 0.0, 2.547366e+04, -0.44668285},
        //     {2.5, 0.0, 0.0, 0.0, 0.0, 2.547366e+04, -0.44668285}, // Same coefficients for the single temperature range
        //     "L 6/94"
        // },
        // O2
        {
            "O2",
            {{"O", 2.0}},
            "NASA7",
            {200.0, 1000.0, 6000.0},
            {3.78245636, -2.99673416e-03, 9.84730201e-06, -9.68129509e-09, 3.24372837e-12, -1063.94356, 3.65767573},
            {3.66096065, 6.56365811e-04, -1.41149627e-07, 2.05797935e-11, -1.29913436e-15, -1215.97718, 3.41536279},
            "RUS 89"
        },
        // O
        {
            "O",
            {{"O", 1.0}},
            "NASA7",
            {200.0, 1000.0, 6000.0},
            {3.1682671, -3.27931884e-03, 6.64306396e-06, -6.12806624e-09, 2.11265971e-12, 2.91222592e+04, 2.05193346},
            {2.54363697, -2.73162486e-05, -4.1902952e-09, 4.95481845e-12, -4.79553694e-16, 2.9226012e+04, 4.92229457},
            "L 1/90"
        },
        // CH2
        {
            "CH2",
            {{"C", 1.0}, {"H", 2.0}},
            "NASA7",
            {200.0, 1000.0, 6000.0},
            {3.71757846, 1.2739126e-03, 2.17347251e-06, -3.488585e-09, 1.65208866e-12, 4.58723866e+04, 1.75297945},
            {3.14631886, 3.03671259e-03, -9.96474439e-07, 1.5048358e-10, -8.57335515e-15, 4.60412605e+04, 4.72341711},
            "IU3/03"
        },
        // CH2(S)
        {
            "CH2(S)",
            {{"C", 1.0}, {"H", 2.0}},
            "NASA7",
            {200.0, 1000.0, 6000.0},
            {4.19331325, -2.33105184e-03, 8.15676451e-06, -6.62985981e-09, 1.93233199e-12, 5.03662246e+04, -0.74673431},
            {3.13501686, 2.89593926e-03, -8.1666809e-07, 1.13572697e-10, -6.36262835e-15, 5.05040504e+04, 4.06030621},
            "IU6/03"
        },
        // N2
        {
            "N2",
            {{"N", 2.0}},
            "NASA7",
            {200.0, 1000.0, 6000.0},
            {3.53100528, -1.23660988e-04, -5.02999433e-07, 2.43530612e-09, -1.40881235e-12, -1046.97628, 2.96747038},
            {2.95257637, 1.3969004e-03, -4.92631603e-07, 7.86010195e-11, -4.60755204e-15, -923.948688, 5.87188762},
            "G 8/02"
        }
    };

    // Test each expected thermodynamic data
    for (const auto& expected : expectedThermo) {
        // Find the corresponding species in the mechanism
        auto it = std::find_if(mechanism.thermoSpecies.begin(), mechanism.thermoSpecies.end(),
            [&expected](const ChemistryVars::ThermoData& thermo) { return thermo.name == expected.name; });

        if (it == mechanism.thermoSpecies.end()) {
            results.totalTests++;
            results.failureMessages.push_back("Species not found: " + expected.name);
            std::cout << " Species not found: " << expected.name << std::endl;
            continue;
        }

        const ChemistryVars::ThermoData& actual = *it;
        bool speciesPassed = true;
        results.totalTests++;

        // Check name
        if (actual.name != expected.name) {
            speciesPassed = false;
            results.failureMessages.push_back("Species name mismatch: " + expected.name);
            std::cout << " Species name mismatch: " << expected.name << std::endl;
        }

        // Check model
        if (actual.model != expected.model) {
            speciesPassed = false;
            results.failureMessages.push_back(expected.name + ": Model mismatch");
            std::cout << " " << expected.name << ": Model mismatch, Expected: " << expected.model
                << ", Actual: " << actual.model << std::endl;
        }

        // Check composition
        for (const auto& element_count : expected.composition) {
            const std::string& element = element_count.first;
            double count = element_count.second;

            if (actual.composition.find(element) == actual.composition.end() ||
                !isEqual(actual.composition.at(element), count)) {
                speciesPassed = false;
                results.failureMessages.push_back(expected.name + ": Composition mismatch");
                std::cout << " " << expected.name << ": Composition mismatch, Element: " << element
                    << ", Expected: " << count << ", Actual: "
                    << (actual.composition.find(element) != actual.composition.end() ?
                        std::to_string(actual.composition.at(element)) : "missing") << std::endl;
            }
        }

        // Check temperature ranges
        if (!compareVectors(actual.temperatureRanges, expected.temperatureRanges)) {
            speciesPassed = false;
            results.failureMessages.push_back(expected.name + ": Temperature range mismatch");
            std::cout << "" << expected.name << ": Temperature range mismatch" << std::endl;
            std::cout << "   Expected: ";
            for (auto t : expected.temperatureRanges) std::cout << t << " ";
            std::cout << std::endl;
            std::cout << "   Actual: ";
            for (auto t : actual.temperatureRanges) std::cout << t << " ";
            std::cout << std::endl;
        }

        // Check coefficients
        if (!compareVectors(actual.coefficients.low, expected.lowCoeffs)) {
            speciesPassed = false;
            results.failureMessages.push_back(expected.name + ": Low temperature coefficients mismatch");
            std::cout << "" << expected.name << ": Low temperature coefficients mismatch" << std::endl;
        }

        if (!compareVectors(actual.coefficients.high, expected.highCoeffs)) {
            speciesPassed = false;
            results.failureMessages.push_back(expected.name + ": High temperature coefficients mismatch");
            std::cout << "" << expected.name << ": High temperature coefficients mismatch" << std::endl;
        }

        if (speciesPassed) {
            results.passedTests++;
            std::cout << " Species " << expected.name << " thermodynamic data test passed" << std::endl;
        }
    }
}

// Test reaction data
void testReactions(ChemistryVars::MechanismData& mechanism, TestResults& results) {
    std::cout << "\n===== Reaction Data Tests =====\n";

    // Expected reaction data
    struct ExpectedReaction {
        std::string equation;
        double A;
        double b;
        double Ea;
        bool isPlogReaction;
        bool isFalloff;
        std::map<std::string, double> efficiencies;
        struct {
            double A;
            double T3;
            double T1;
            double T2;
        } troe;
        struct {
            double A;
            double b;
            double Ea;
        } lowRate;
    };

    std::vector<ExpectedReaction> expectedReactions = {
        // CH3 + O2 <=> CH2O + OH
        {
            "CH3 + O2 <=> CH2O + OH",
            2.641,
            3.283,
            8105.0,
            false,
            false,
            {},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0}
        },
        // CH3 + O <=> CH2O + H
        {
            "CH3 + O <=> CH2O + H",
            5.54e+13,
            0.05,
            -136.0,
            false,
            false,
            {},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0}
        },
        // CH3 + OH <=> CH2(S) + H2O
        {
            "CH3 + OH <=> CH2(S) + H2O",
            5.282e+17,
            -1.518,
            1772.0,
            false,
            false,
            {},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0}
        },
        // CH3 + OH <=> CH2O + H2 (PLOG)
        {
            "CH3 + OH <=> CH2O + H2",
            1.65e+07,
            0.973,
            -2010.0,
            true,
            false,
            {},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0}
        },
        // CH3OH (+M) <=> CH3 + OH (+M) (falloff)
        {
            "CH3OH (+M) <=> CH3 + OH (+M)",
            2.084e+18,
            -0.615,
            9.25406e+04,
            false,
            true,
            {},
            {-0.4748, 3.558e+04, 1116.0, 9023.0},
            {1.5e+43, -6.995, 9.79922e+04}
        }
    };

    // Create a copy of reactions list so we can remove matched reactions
    std::vector<ChemistryVars::ReactionData> availableReactions = mechanism.reactions;

    // Test each expected reaction
    for (const auto& expected : expectedReactions) {
        // Find all reactions with matching equation
        std::vector<size_t> matchingIndices;
        for (size_t i = 0; i < availableReactions.size(); i++) {
            if (availableReactions[i].equation == expected.equation) {
                matchingIndices.push_back(i);
            }
        }

        if (matchingIndices.empty()) {
            results.totalTests++;
            results.failureMessages.push_back("Reaction not found: " + expected.equation);
            std::cout << "Reaction not found: " << expected.equation << std::endl;
            continue;
        }

        // If we found multiple matching reactions, find the best match by rate constants
        size_t bestMatchIndex = matchingIndices[0];
        double bestMatchScore = std::numeric_limits<double>::max();

        for (size_t idx : matchingIndices) {
            const auto& rxn = availableReactions[idx];

            // Calculate a score based on differences (lower is better match)
            double scoreA = std::abs(rxn.rateConstant.A - expected.A) / (expected.A > 1.0 ? expected.A : 1.0);
            double scoreB = std::abs(rxn.rateConstant.b - expected.b);
            double scoreEa = std::abs(rxn.rateConstant.Ea - expected.Ea) /
                (std::abs(expected.Ea) > 1.0 ? std::abs(expected.Ea) : 1.0);

            double score = scoreA + scoreB + scoreEa;

            if (score < bestMatchScore) {
                bestMatchScore = score;
                bestMatchIndex = idx;
            }
        }

        // Use the best match for testing
        const ChemistryVars::ReactionData& actual = availableReactions[bestMatchIndex];
        bool reactionPassed = true;
        results.totalTests++;

        // Check rate constants
        double Atol = expected.A * 1e-4 > 1e-10 ? expected.A * 1e-4 : 1e-10;
        double Eatol = std::abs(expected.Ea) * 1e-4 > 1.0 ? std::abs(expected.Ea) * 1e-4 : 1.0;

        if (!isEqual(actual.rateConstant.A, expected.A, Atol) ||
            !isEqual(actual.rateConstant.b, expected.b, 0.05) ||
            !isEqual(actual.rateConstant.Ea, expected.Ea, Eatol)) {
            reactionPassed = false;
            results.failureMessages.push_back(expected.equation + ": Rate constants mismatch");
            std::cout << " " << expected.equation << ": Rate constants mismatch" << std::endl;
            std::cout << "   Expected: A=" << expected.A << ", b=" << expected.b << ", Ea=" << expected.Ea << std::endl;
            std::cout << "   Actual: A=" << actual.rateConstant.A << ", b=" << actual.rateConstant.b
                << ", Ea=" << actual.rateConstant.Ea << std::endl;
        }

        // Check PLOG flag - only when expected is a PLOG reaction
        if (expected.isPlogReaction && !actual.rateConstant.isPressureDependent) {
            reactionPassed = false;
            results.failureMessages.push_back(expected.equation + ": Missing PLOG flag");
            std::cout << "" << expected.equation << ": Should be a PLOG reaction but not marked as pressure dependent" << std::endl;
        }

        // Check falloff reactions with low pressure data and Troe parameters
        if (expected.isFalloff) {
            std::string actualType = actual.type;
            // Check if reaction is of Falloff type
            if (actualType != "falloff" && actualType != "chemically-activated" && actualType != "three-body") {
                reactionPassed = false;
                results.failureMessages.push_back(expected.equation + ": Missing Falloff flag");
                std::cout << " " << expected.equation << ": Should be a Falloff reaction but is " << actualType << std::endl;
                continue;
            }

            // Check low-pressure data
            double lowAtol = expected.lowRate.A * 1e-4 > 1e-10 ? expected.lowRate.A * 1e-4 : 1e-10;
            double lowEatol = std::abs(expected.lowRate.Ea) * 1e-4 > 1.0 ? std::abs(expected.lowRate.Ea) * 1e-4 : 1.0;

            if (!isEqual(actual.lowPressure.A, expected.lowRate.A, lowAtol) ||
                !isEqual(actual.lowPressure.b, expected.lowRate.b, 0.05) ||
                !isEqual(actual.lowPressure.Ea, expected.lowRate.Ea, lowEatol)) {
                reactionPassed = false;
                results.failureMessages.push_back(expected.equation + ": Low pressure data mismatch");
                std::cout << " " << expected.equation << ": Low pressure data mismatch" << std::endl;
                std::cout << "   Expected: A=" << expected.lowRate.A << ", b=" << expected.lowRate.b
                    << ", Ea=" << expected.lowRate.Ea << std::endl;
                std::cout << "   Actual: A=" << actual.lowPressure.A << ", b=" << actual.lowPressure.b
                    << ", Ea=" << actual.lowPressure.Ea << std::endl;
            }

            // Check Troe parameters (if present)
            if (expected.troe.A != 0.0 || expected.troe.T1 != 0.0 || expected.troe.T2 != 0.0 || expected.troe.T3 != 0.0) {
                if (!isEqual(actual.troe.a, expected.troe.A, 0.01) ||
                    !isEqual(actual.troe.T_star, expected.troe.T3, expected.troe.T3 * 0.01) ||
                    !isEqual(actual.troe.T_double_star, expected.troe.T1, expected.troe.T1 * 0.01) ||
                    !isEqual(actual.troe.T_triple_star, expected.troe.T2, expected.troe.T2 * 0.01)) {
                    reactionPassed = false;
                    results.failureMessages.push_back(expected.equation + ": Troe parameters mismatch");
                    std::cout << " " << expected.equation << ": Troe parameters mismatch" << std::endl;
                    std::cout << "   Expected: A=" << expected.troe.A << ", T3=" << expected.troe.T3
                        << ", T1=" << expected.troe.T1 << ", T2=" << expected.troe.T2 << std::endl;
                    std::cout << "   Actual: A=" << actual.troe.a << ", T3=" << actual.troe.T_star
                        << ", T1=" << actual.troe.T_double_star << ", T2=" << actual.troe.T_triple_star << std::endl;
                }
            }

            // Check third-body efficiencies (if present)
            for (const auto& efficiency_pair : expected.efficiencies) {
                const std::string& species = efficiency_pair.first;
                double eff = efficiency_pair.second;

                if (actual.efficiencies.find(species) == actual.efficiencies.end() ||
                    !isEqual(actual.efficiencies.at(species), eff, 0.01)) {
                    reactionPassed = false;
                    results.failureMessages.push_back(expected.equation + ": Third-body efficiency mismatch");
                    std::cout << " " << expected.equation << ": Species " << species << " efficiency mismatch, Expected: "
                        << eff << ", Actual: "
                        << (actual.efficiencies.find(species) != actual.efficiencies.end() ?
                            std::to_string(actual.efficiencies.at(species)) : "missing") << std::endl;
                }
            }
        }

        if (reactionPassed) {
            results.passedTests++;
            std::cout << " Reaction " << expected.equation << " test passed" << std::endl;
        }

        // Remove the matched reaction from available list to avoid matching it again
        availableReactions.erase(availableReactions.begin() + bestMatchIndex);
    }
}

// Test transport data
void testTransport(ChemistryVars::MechanismData& mechanism, TestResults& results) {
    std::cout << "\n===== Transport Data Tests =====\n";

    // Expected transport data
    struct ExpectedTransport {
        std::string name;
        std::string model;
        std::string geometry;
        double diameter;
        double wellDepth;
        double polarizability;
        double rotationalRelaxation;
        std::string note;
    };

    std::vector<ExpectedTransport> expectedTransport = {
        // H2
        {
            "H2",
            "gas",
            "linear",
            2.92,
            38.0,
            0.79,
            280.0,
            "!"
        },
        // O2
        {
            "O2",
            "gas",
            "linear",
            3.458,
            107.4,
            1.6,
            3.8,
            "!"
        },
        // O
        {
            "O",
            "gas",
            "atom",
            2.75,
            80.0,
            0.0, // Default value as original data doesn't provide it
            0.0, // Default value as original data doesn't provide it
            "!"
        },
        // N2
        {
            "N2",
            "gas",
            "linear",
            3.621,
            97.53,
            1.76,
            4.0,
            "!"
        }
    };

    // Test each expected transport data
    for (const auto& expected : expectedTransport) {
        // Find the corresponding species in the mechanism
        auto it = std::find_if(mechanism.transportSpecies.begin(), mechanism.transportSpecies.end(),
            [&expected](const ChemistryVars::TransportData& transport) { return transport.name == expected.name; });

        if (it == mechanism.transportSpecies.end()) {
            results.totalTests++;
            results.failureMessages.push_back("Transport data not found for species: " + expected.name);
            std::cout << " Transport data not found for species: " << expected.name << std::endl;
            continue;
        }

        const ChemistryVars::TransportData& actual = *it;
        bool transportPassed = true;
        results.totalTests++;

        // Check model
        if (actual.model != expected.model) {
            transportPassed = false;
            results.failureMessages.push_back(expected.name + ": Transport model mismatch");
            std::cout << " " << expected.name << ": Transport model mismatch, Expected: " << expected.model
                << ", Actual: " << actual.model << std::endl;
        }

        // Check geometry
        if (actual.geometry != expected.geometry) {
            transportPassed = false;
            results.failureMessages.push_back(expected.name + ": Geometry mismatch");
            std::cout << " " << expected.name << ": Geometry mismatch, Expected: " << expected.geometry
                << ", Actual: " << actual.geometry << std::endl;
        }

        // Check diameter
        if (!isEqual(actual.diameter, expected.diameter, 0.01)) {
            transportPassed = false;
            results.failureMessages.push_back(expected.name + ": Diameter mismatch");
            std::cout << " " << expected.name << ": Diameter mismatch, Expected: " << expected.diameter
                << ", Actual: " << actual.diameter << std::endl;
        }

        // Check well depth
        if (!isEqual(actual.wellDepth, expected.wellDepth, 0.1)) {
            transportPassed = false;
            results.failureMessages.push_back(expected.name + ": Well depth mismatch");
            std::cout << " " << expected.name << ": Well depth mismatch, Expected: " << expected.wellDepth
                << ", Actual: " << actual.wellDepth << std::endl;
        }

        // Check polarizability (only if expected value is non-zero)
        if (expected.polarizability > 0.0 && !isEqual(actual.polarizability, expected.polarizability, 0.01)) {
            transportPassed = false;
            results.failureMessages.push_back(expected.name + ": Polarizability mismatch");
            std::cout << " " << expected.name << ": Polarizability mismatch, Expected: " << expected.polarizability
                << ", Actual: " << actual.polarizability << std::endl;
        }

        // Check rotational relaxation (only if expected value is non-zero)
        if (expected.rotationalRelaxation > 0.0 && !isEqual(actual.rotationalRelaxation, expected.rotationalRelaxation, 0.1)) {
            transportPassed = false;
            results.failureMessages.push_back(expected.name + ": Rotational relaxation mismatch");
            std::cout << " " << expected.name << ": Rotational relaxation mismatch, Expected: " << expected.rotationalRelaxation
                << ", Actual: " << actual.rotationalRelaxation << std::endl;
        }

        if (transportPassed) {
            results.passedTests++;
            std::cout << " Species " << expected.name << " transport data test passed" << std::endl;
        }
    }
}

// Main test function
bool runMechanismTest(const std::string& yamlFile) {
    std::cout << "Starting general test using file: " << yamlFile << std::endl;

    // Initialize test results
    TestResults results;

    try {
        // Load mechanism data
        std::cout << "Loading mechanism data..." << std::endl;
        ChemistryVars::MechanismData mechanism = ChemistryVars::loadMechanism(yamlFile, false);

        std::cout << "Data loaded successfully, containing:" << std::endl
            << " - " << mechanism.reactions.size() << " reactions" << std::endl
            << " - " << mechanism.thermoSpecies.size() << " thermodynamic species" << std::endl
            << " - " << mechanism.transportSpecies.size() << " transport data entries" << std::endl;

        // Run tests
        testThermo(mechanism, results);
        testReactions(mechanism, results);
        testTransport(mechanism, results);

        // Print test result summary
        std::cout << "\n===== Test Results Summary =====\n";
        std::cout << "Total tests: " << results.totalTests << std::endl;
        std::cout << "Passed tests: " << results.passedTests << std::endl;
        std::cout << "Success rate: " << std::fixed << std::setprecision(2)
            << (results.totalTests > 0 ? 100.0 * results.passedTests / results.totalTests : 0)
            << "%" << std::endl;

        if (results.failureMessages.size() > 0) {
            std::cout << "\nFailed tests:\n";
            for (const auto& msg : results.failureMessages) {
                std::cout << " - " << msg << std::endl;
            }
        }

        return results.passedTests == results.totalTests;
    }
    catch (const std::exception& e) {
        std::cerr << "Error during test execution: " << e.what() << std::endl;
        return false;
    }
}

// 测试特定YAML数据
bool testSpecificC2H4Data(const std::string& yamlFile) {
    std::cout << "Starting C2H4 specific data test using file: " << yamlFile << std::endl;

    // Initialize test results
    TestResults results;

    try {
        // Load mechanism data
        std::cout << "Loading mechanism data..." << std::endl;
        ChemistryVars::MechanismData mechanism = ChemistryVars::loadMechanism(yamlFile, false);

        std::cout << "Data loaded successfully, starting specific data tests..." << std::endl;

        // 1. Test specific species thermodynamic data
        std::cout << "\n===== Specific Species Thermodynamic Data Tests =====\n";

        // First test H species
        {
            const std::string speciesName = "H";
            auto it = std::find_if(mechanism.thermoSpecies.begin(), mechanism.thermoSpecies.end(),
                [&speciesName](const ChemistryVars::ThermoData& thermo) { return thermo.name == speciesName; });

            if (it == mechanism.thermoSpecies.end()) {
                results.totalTests++;
                results.failureMessages.push_back("Species not found: " + speciesName);
                std::cout << " Species not found: " << speciesName << std::endl;
            }
            else {
                const ChemistryVars::ThermoData& actual = *it;
                bool passed = true;
                results.totalTests++;

                // Verify model
                if (actual.model != "NASA7") {
                    passed = false;
                    results.failureMessages.push_back(speciesName + ": Model mismatch");
                    std::cout << " " << speciesName << ": Model mismatch, Expected: NASA7, Actual: " << actual.model << std::endl;
                }

                // Verify temperature range (H has only one range here)
                std::vector<double> expectedTempRanges = { 200.0, 6000.0 };
                if (actual.temperatureRanges.size() != 2 ||
                    !isEqual(actual.temperatureRanges[0], 200.0) ||
                    !isEqual(actual.temperatureRanges[1], 6000.0)) {
                    passed = false;
                    results.failureMessages.push_back(speciesName + ": Temperature range mismatch");
                    std::cout << " " << speciesName << ": Temperature range mismatch" << std::endl;
                    std::cout << "   Expected: 200.0 6000.0" << std::endl;
                    std::cout << "   Actual: ";
                    for (auto t : actual.temperatureRanges) std::cout << t << " ";
                    std::cout << std::endl;
                }

                // Verify composition
                if (actual.composition.size() != 1 ||
                    actual.composition.find("H") == actual.composition.end() ||
                    !isEqual(actual.composition.at("H"), 1.0)) {
                    passed = false;
                    results.failureMessages.push_back(speciesName + ": Composition mismatch");
                    std::cout << " " << speciesName << ": Composition mismatch, Expected: H:1, Actual: ";
                    //for (const auto& [element, count] : actual.composition) {
                    //    std::cout << element << ":" << count << " ";
                    //}
                    for (const auto& pair : actual.composition) {
                        const auto& element = pair.first;
                        const auto& count = pair.second;
                        std::cout << element << ":" << count << " ";
                    }
                    std::cout << std::endl;
                }

                // Verify coefficients (H has only one set of coefficients)
                std::vector<double> expectedCoeffs = { 2.5, 0.0, 0.0, 0.0, 0.0, 2.547366e+04, -0.44668285 };
                if (!compareVectors(actual.coefficients.low, expectedCoeffs)) {
                    passed = false;
                    results.failureMessages.push_back(speciesName + ": Coefficients mismatch");
                    std::cout << " " << speciesName << ": Coefficients mismatch" << std::endl;
                    std::cout << "   Expected: 2.5 0.0 0.0 0.0 0.0 2.547366e+04 -0.44668285" << std::endl;
                    std::cout << "   Actual: ";
                    for (auto c : actual.coefficients.low) std::cout << c << " ";
                    std::cout << std::endl;
                }

                if (passed) {
                    results.passedTests++;
                    std::cout << " Species " << speciesName << " thermodynamic data test passed" << std::endl;
                }
            }
        }

        // Test HCO species
        {
            const std::string speciesName = "HCO";
            auto it = std::find_if(mechanism.thermoSpecies.begin(), mechanism.thermoSpecies.end(),
                [&speciesName](const ChemistryVars::ThermoData& thermo) { return thermo.name == speciesName; });

            if (it == mechanism.thermoSpecies.end()) {
                results.totalTests++;
                results.failureMessages.push_back("Species not found: " + speciesName);
                std::cout << " Species not found: " << speciesName << std::endl;
            }
            else {
                const ChemistryVars::ThermoData& actual = *it;
                bool passed = true;
                results.totalTests++;

                // Verify composition
                std::map<std::string, double> expectedComp = { {"C", 1.0}, {"H", 1.0}, {"O", 1.0} };
                //for (const auto& [element, count] : expectedComp) 
                //{
                for (const auto& pair : expectedComp) {
                    const auto& element = pair.first;
                    const auto& count = pair.second;
                    if (actual.composition.find(element) == actual.composition.end() ||
                        !isEqual(actual.composition.at(element), count)) {
                        passed = false;
                        results.failureMessages.push_back(speciesName + ": Composition mismatch");
                        std::cout << " " << speciesName << ": Element " << element << " mismatch, Expected: " << count;
                        if (actual.composition.find(element) != actual.composition.end()) {
                            std::cout << ", Actual: " << actual.composition.at(element);
                        }
                        else {
                            std::cout << ", Actual: missing";
                        }
                        std::cout << std::endl;
                        break;
                    }
                }

                // Verify temperature range
                std::vector<double> expectedTempRanges = { 200.0, 1000.0, 6000.0 };
                if (!compareVectors(actual.temperatureRanges, expectedTempRanges)) {
                    passed = false;
                    results.failureMessages.push_back(speciesName + ": Temperature range mismatch");
                    std::cout << " " << speciesName << ": Temperature range mismatch" << std::endl;
                    std::cout << "   Expected: 200.0 1000.0 6000.0" << std::endl;
                    std::cout << "   Actual: ";
                    for (auto t : actual.temperatureRanges) std::cout << t << " ";
                    std::cout << std::endl;
                }

                // Verify low temperature coefficients
                std::vector<double> expectedLowCoeffs = { 4.2375461, -3.32075257e-03, 1.40030264e-05, -1.34239995e-08, 4.37416208e-12, 3872.41185, 3.30834869 };
                if (!compareVectors(actual.coefficients.low, expectedLowCoeffs, 1e-4)) {
                    passed = false;
                    results.failureMessages.push_back(speciesName + ": Low temperature coefficients mismatch");
                    std::cout << " " << speciesName << ": Low temperature coefficients mismatch" << std::endl;
                }

                // Verify high temperature coefficients
                std::vector<double> expectedHighCoeffs = { 3.92001542, 2.52279324e-03, -6.71004164e-07, 1.05615948e-10, -7.43798261e-15, 3653.42928, 3.58077056 };
                if (!compareVectors(actual.coefficients.high, expectedHighCoeffs, 1e-4)) {
                    passed = false;
                    results.failureMessages.push_back(speciesName + ": High temperature coefficients mismatch");
                    std::cout << " " << speciesName << ": High temperature coefficients mismatch" << std::endl;
                }

                if (passed) {
                    results.passedTests++;
                    std::cout << " Species " << speciesName << " thermodynamic data test passed" << std::endl;
                }
            }
        }

        // Test HO2C2H3OH species (single temperature range with lower coefficients)
        {
            const std::string speciesName = "HO2C2H3OH";
            auto it = std::find_if(mechanism.thermoSpecies.begin(), mechanism.thermoSpecies.end(),
                [&speciesName](const ChemistryVars::ThermoData& thermo) { return thermo.name == speciesName; });

            if (it == mechanism.thermoSpecies.end()) {
                results.totalTests++;
                results.failureMessages.push_back("Species not found: " + speciesName);
                std::cout << " Species not found: " << speciesName << std::endl;
            }
            else {
                const ChemistryVars::ThermoData& actual = *it;
                bool passed = true;
                results.totalTests++;

                // Verify temperature range (special case)
                std::vector<double> expectedTempRanges = { 300.0, 2000.0 };
                if (!compareVectors(actual.temperatureRanges, expectedTempRanges)) {
                    passed = false;
                    results.failureMessages.push_back(speciesName + ": Temperature range mismatch");
                    std::cout << " " << speciesName << ": Temperature range mismatch" << std::endl;
                    std::cout << "   Expected: 300.0 2000.0" << std::endl;
                    std::cout << "   Actual: ";
                    for (auto t : actual.temperatureRanges) std::cout << t << " ";
                    std::cout << std::endl;
                }

                // Verify composition
                std::map<std::string, double> expectedComp = { {"C", 2.0}, {"H", 5.0}, {"O", 3.0} };
                //for (const auto& [element, count] : expectedComp) {
                for (const auto& pair : expectedComp) {
                    const auto& element = pair.first;
                    const auto& count = pair.second;
                    if (actual.composition.find(element) == actual.composition.end() ||
                        !isEqual(actual.composition.at(element), count)) {
                        passed = false;
                        results.failureMessages.push_back(speciesName + ": Composition mismatch");
                        std::cout << " " << speciesName << ": Element " << element << " mismatch, Expected: " << count;
                        if (actual.composition.find(element) != actual.composition.end()) {
                            std::cout << ", Actual: " << actual.composition.at(element);
                        }
                        else {
                            std::cout << ", Actual: missing";
                        }
                        std::cout << std::endl;
                        break;
                    }
                }

                // Verify has low coefficients but no high coefficients (special case)
                if (actual.coefficients.low.empty()) {
                    passed = false;
                    results.failureMessages.push_back(speciesName + ": Missing coefficients");
                    std::cout << " " << speciesName << ": Missing coefficients" << std::endl;
                }

                if (!actual.coefficients.high.empty()) {
                    passed = false;
                    results.failureMessages.push_back(speciesName + ": Should not have high coefficients");
                    std::cout << " " << speciesName << ": Should not have high coefficients" << std::endl;
                }

                // Verify coefficient values
                std::vector<double> expectedCoeffs = { 17.6346697, -0.071859728, 2.55608152e-04, -3.01649715e-07, 1.21090555e-10, -2.02049117e+04, -43.4862383 };
                if (!compareVectors(actual.coefficients.low, expectedCoeffs, 1e-4)) {
                    passed = false;
                    results.failureMessages.push_back(speciesName + ": Coefficients mismatch");
                    std::cout << " " << speciesName << ": Coefficients mismatch" << std::endl;
                }

                if (passed) {
                    results.passedTests++;
                    std::cout << " Species " << speciesName << " thermodynamic data test passed" << std::endl;
                }
            }
        }

        // 2. Test specific reaction data
        std::cout << "\n===== Specific Reaction Data Tests =====\n";

        // Test simple reaction: O + H2O <=> OH + OH
        {
            const std::string rxnEquation = "O + H2O <=> OH + OH";
            auto it = std::find_if(mechanism.reactions.begin(), mechanism.reactions.end(),
                [&rxnEquation](const ChemistryVars::ReactionData& rxn) { return rxn.equation == rxnEquation; });

            if (it == mechanism.reactions.end()) {
                results.totalTests++;
                results.failureMessages.push_back("Reaction not found: " + rxnEquation);
                std::cout << " Reaction not found: " << rxnEquation << std::endl;
            }
            else {
                const ChemistryVars::ReactionData& actual = *it;
                bool passed = true;
                results.totalTests++;

                // Verify rate constants
                if (!isEqual(actual.rateConstant.A, 6.7e+07, 6.7e+05) ||
                    !isEqual(actual.rateConstant.b, 1.704, 0.01) ||
                    !isEqual(actual.rateConstant.Ea, 1.49868e+04, 10.0)) {
                    passed = false;
                    results.failureMessages.push_back(rxnEquation + ": Rate constants mismatch");
                    std::cout << " " << rxnEquation << ": Rate constants mismatch" << std::endl;
                    std::cout << "   Expected: A=6.7e+07, b=1.704, Ea=1.49868e+04" << std::endl;
                    std::cout << "   Actual: A=" << actual.rateConstant.A << ", b=" << actual.rateConstant.b
                        << ", Ea=" << actual.rateConstant.Ea << std::endl;
                }

                if (passed) {
                    results.passedTests++;
                    std::cout << " Reaction " << rxnEquation << " test passed" << std::endl;
                }
            }
        }

        // Test three-body reaction: O + H + M <=> OH + M
        {
            const std::string rxnEquation = "O + H + M <=> OH + M";
            auto it = std::find_if(mechanism.reactions.begin(), mechanism.reactions.end(),
                [&rxnEquation](const ChemistryVars::ReactionData& rxn) { return rxn.equation == rxnEquation; });

            if (it == mechanism.reactions.end()) {
                results.totalTests++;
                results.failureMessages.push_back("Reaction not found: " + rxnEquation);
                std::cout << " Reaction not found: " << rxnEquation << std::endl;
            }
            else {
                const ChemistryVars::ReactionData& actual = *it;
                bool passed = true;
                results.totalTests++;

                // Verify reaction type
                if (actual.type != "three-body") {
                    passed = false;
                    results.failureMessages.push_back(rxnEquation + ": Reaction type mismatch");
                    std::cout << " " << rxnEquation << ": Reaction type mismatch, Expected: three-body, Actual: " << actual.type << std::endl;
                }

                // Verify rate constants
                if (!isEqual(actual.rateConstant.A, 4.714e+18, 4.714e+16) ||
                    !isEqual(actual.rateConstant.b, -1.0, 0.01) ||
                    !isEqual(actual.rateConstant.Ea, 0.0, 0.1)) {
                    passed = false;
                    results.failureMessages.push_back(rxnEquation + ": Rate constants mismatch");
                    std::cout << " " << rxnEquation << ": Rate constants mismatch" << std::endl;
                    std::cout << "   Expected: A=4.714e+18, b=-1.0, Ea=0.0" << std::endl;
                    std::cout << "   Actual: A=" << actual.rateConstant.A << ", b=" << actual.rateConstant.b
                        << ", Ea=" << actual.rateConstant.Ea << std::endl;
                }

                // Verify third-body efficiencies
                std::map<std::string, double> expectedEfficiencies = {
                    {"H2", 2.5}, {"H2O", 12.0}, {"CO", 1.5}, {"CO2", 2.0}, {"CH4", 2.0}
                };

                //for (const auto& [species, eff] : expectedEfficiencies) {
                for (const auto& pair : expectedEfficiencies) {
                    const auto& species = pair.first;
                    const auto& eff = pair.second;
                    if (actual.efficiencies.find(species) == actual.efficiencies.end() ||
                        !isEqual(actual.efficiencies.at(species), eff, 0.01)) {
                        passed = false;
                        results.failureMessages.push_back(rxnEquation + ": Third-body efficiency mismatch");
                        std::cout << " " << rxnEquation << ": Species " << species << " efficiency mismatch, Expected: " << eff;
                        if (actual.efficiencies.find(species) != actual.efficiencies.end()) {
                            std::cout << ", Actual: " << actual.efficiencies.at(species);
                        }
                        else {
                            std::cout << ", Actual: missing";
                        }
                        std::cout << std::endl;
                        break;
                    }
                }

                if (passed) {
                    results.passedTests++;
                    std::cout << " Reaction " << rxnEquation << " test passed" << std::endl;
                }
            }
        }

        // Test Falloff reaction: H2O2 (+M) <=> OH + OH (+M)
        {
            const std::string rxnEquation = "H2O2 (+M) <=> OH + OH (+M)";
            auto it = std::find_if(mechanism.reactions.begin(), mechanism.reactions.end(),
                [&rxnEquation](const ChemistryVars::ReactionData& rxn) { return rxn.equation == rxnEquation; });

            if (it == mechanism.reactions.end()) {
                results.totalTests++;
                results.failureMessages.push_back("Reaction not found: " + rxnEquation);
                std::cout << " Reaction not found: " << rxnEquation << std::endl;
            }
            else {
                const ChemistryVars::ReactionData& actual = *it;
                bool passed = true;
                results.totalTests++;

                // Verify reaction type
                if (actual.type != "falloff") {
                    passed = false;
                    results.failureMessages.push_back(rxnEquation + ": Reaction type mismatch");
                    std::cout << " " << rxnEquation << ": Reaction type mismatch, Expected: falloff, Actual: " << actual.type << std::endl;
                }

                // Verify high-pressure rate constants
                if (!isEqual(actual.rateConstant.A, 2.0e+12, 2.0e+10) ||
                    !isEqual(actual.rateConstant.b, 0.9, 0.01) ||
                    !isEqual(actual.rateConstant.Ea, 4.8749e+04, 50.0)) {
                    passed = false;
                    results.failureMessages.push_back(rxnEquation + ": High-pressure rate constants mismatch");
                    std::cout << " " << rxnEquation << ": High-pressure rate constants mismatch" << std::endl;
                    std::cout << "   Expected: A=2.0e+12, b=0.9, Ea=4.8749e+04" << std::endl;
                    std::cout << "   Actual: A=" << actual.rateConstant.A << ", b=" << actual.rateConstant.b
                        << ", Ea=" << actual.rateConstant.Ea << std::endl;
                }

                // Verify low-pressure rate constants
                if (!isEqual(actual.lowPressure.A, 2.49e+24, 2.49e+22) ||
                    !isEqual(actual.lowPressure.b, -2.3, 0.01) ||
                    !isEqual(actual.lowPressure.Ea, 4.8749e+04, 50.0)) {
                    passed = false;
                    results.failureMessages.push_back(rxnEquation + ": Low-pressure rate constants mismatch");
                    std::cout << " " << rxnEquation << ": Low-pressure rate constants mismatch" << std::endl;
                    std::cout << "   Expected: A=2.49e+24, b=-2.3, Ea=4.8749e+04" << std::endl;
                    std::cout << "   Actual: A=" << actual.lowPressure.A << ", b=" << actual.lowPressure.b
                        << ", Ea=" << actual.lowPressure.Ea << std::endl;
                }

                // Verify Troe parameters
                if (!isEqual(actual.troe.a, 0.43, 0.01) ||
                    !isEqual(actual.troe.T_star, 1.0e-20, 1.0e-22) ||
                    !isEqual(actual.troe.T_double_star, 1.0e+20, 1.0e+18)) {
                    passed = false;
                    results.failureMessages.push_back(rxnEquation + ": Troe parameters mismatch");
                    std::cout << " " << rxnEquation << ": Troe parameters mismatch" << std::endl;
                    std::cout << "   Expected: a=0.43, T3=1.0e-20, T1=1.0e+20" << std::endl;
                    std::cout << "   Actual: a=" << actual.troe.a << ", T3=" << actual.troe.T_star
                        << ", T1=" << actual.troe.T_double_star << std::endl;
                }

                // Verify third-body efficiencies
                std::map<std::string, double> expectedEfficiencies = {
                    {"H2O", 7.65}, {"H2", 3.7}, {"CO", 2.8}
                };

                //for (const auto& [species, eff] : expectedEfficiencies) {
                for (const auto& pair : expectedEfficiencies) {
                    const auto& species = pair.first;
                    const auto& eff = pair.second;
                    if (actual.efficiencies.find(species) == actual.efficiencies.end() ||
                        !isEqual(actual.efficiencies.at(species), eff, 0.01)) {
                        passed = false;
                        results.failureMessages.push_back(rxnEquation + ": Third-body efficiency mismatch");
                        std::cout << " " << rxnEquation << ": Species " << species << " efficiency mismatch, Expected: " << eff;
                        if (actual.efficiencies.find(species) != actual.efficiencies.end()) {
                            std::cout << ", Actual: " << actual.efficiencies.at(species);
                        }
                        else {
                            std::cout << ", Actual: missing";
                        }
                        std::cout << std::endl;
                        break;
                    }
                }

                if (passed) {
                    results.passedTests++;
                    std::cout << " Reaction " << rxnEquation << " test passed" << std::endl;
                }
            }
        }

        // Test PLOG reaction: CH3 + OH <=> CH2 + H2O
        {
            const std::string rxnEquation = "CH3 + OH <=> CH2 + H2O";
            auto it = std::find_if(mechanism.reactions.begin(), mechanism.reactions.end(),
                [&rxnEquation](const ChemistryVars::ReactionData& rxn) { return rxn.equation == rxnEquation; });

            if (it == mechanism.reactions.end()) {
                results.totalTests++;
                results.failureMessages.push_back("Reaction not found: " + rxnEquation);
                std::cout << " Reaction not found: " << rxnEquation << std::endl;
            }
            else {
                const ChemistryVars::ReactionData& actual = *it;
                bool passed = true;
                results.totalTests++;

                // Verify rate constants
                if (!isEqual(actual.rateConstant.A, 4.293e+04, 4.293e+02) ||
                    !isEqual(actual.rateConstant.b, 2.568, 0.01) ||
                    !isEqual(actual.rateConstant.Ea, 3997.8, 10.0)) {
                    passed = false;
                    results.failureMessages.push_back(rxnEquation + ": Rate constants mismatch");
                    std::cout << " " << rxnEquation << ": Rate constants mismatch" << std::endl;
                    std::cout << "   Expected: A=4.293e+04, b=2.568, Ea=3997.8" << std::endl;
                    std::cout << "   Actual: A=" << actual.rateConstant.A << ", b=" << actual.rateConstant.b
                        << ", Ea=" << actual.rateConstant.Ea << std::endl;
                }

                // Verify if it's a PLOG reaction
                if (!actual.rateConstant.isPressureDependent) {
                    passed = false;
                    results.failureMessages.push_back(rxnEquation + ": Should be a PLOG reaction");
                    std::cout << " " << rxnEquation << ": Should be a PLOG reaction but not marked as pressure dependent" << std::endl;
                }

                // Verify PLOG data count
                if (actual.rateConstant.plogData.size() < 3) {
                    passed = false;
                    results.failureMessages.push_back(rxnEquation + ": Insufficient PLOG data");
                    std::cout << " " << rxnEquation << ": Insufficient PLOG data, Expected at least 3 sets, Actual: "
                        << actual.rateConstant.plogData.size() << std::endl;
                }

                if (passed) {
                    results.passedTests++;
                    std::cout << " Reaction " << rxnEquation << " test passed" << std::endl;
                }
            }
        }

        // Test C2H3 + H (+M) <=> C2H4 (+M) reaction
        {
            const std::string rxnEquation = "C2H3 + H (+M) <=> C2H4 (+M)";
            auto it = std::find_if(mechanism.reactions.begin(), mechanism.reactions.end(),
                [&rxnEquation](const ChemistryVars::ReactionData& rxn) { return rxn.equation == rxnEquation; });

            if (it == mechanism.reactions.end()) {
                results.totalTests++;
                results.failureMessages.push_back("Reaction not found: " + rxnEquation);
                std::cout << " Reaction not found: " << rxnEquation << std::endl;
            }
            else {
                const ChemistryVars::ReactionData& actual = *it;
                bool passed = true;
                results.totalTests++;

                // Verify reaction type
                if (actual.type != "falloff") {
                    passed = false;
                    results.failureMessages.push_back(rxnEquation + ": Reaction type mismatch");
                    std::cout << " " << rxnEquation << ": Reaction type mismatch, Expected: falloff, Actual: " << actual.type << std::endl;
                }

                // Verify high-pressure rate constants
                if (!isEqual(actual.rateConstant.A, 6.08e+12, 6.08e+10) ||
                    !isEqual(actual.rateConstant.b, 0.27, 0.01) ||
                    !isEqual(actual.rateConstant.Ea, 280.0, 5.0)) {
                    passed = false;
                    results.failureMessages.push_back(rxnEquation + ": High-pressure rate constants mismatch");
                    std::cout << " " << rxnEquation << ": High-pressure rate constants mismatch" << std::endl;
                    std::cout << "   Expected: A=6.08e+12, b=0.27, Ea=280.0" << std::endl;
                    std::cout << "   Actual: A=" << actual.rateConstant.A << ", b=" << actual.rateConstant.b
                        << ", Ea=" << actual.rateConstant.Ea << std::endl;
                }

                // Verify low-pressure rate constants
                if (!isEqual(actual.lowPressure.A, 1.4e+30, 1.4e+28) ||
                    !isEqual(actual.lowPressure.b, -3.86, 0.01) ||
                    !isEqual(actual.lowPressure.Ea, 3320.0, 10.0)) {
                    passed = false;
                    results.failureMessages.push_back(rxnEquation + ": Low-pressure rate constants mismatch");
                    std::cout << " " << rxnEquation << ": Low-pressure rate constants mismatch" << std::endl;
                    std::cout << "   Expected: A=1.4e+30, b=-3.86, Ea=3320.0" << std::endl;
                    std::cout << "   Actual: A=" << actual.lowPressure.A << ", b=" << actual.lowPressure.b
                        << ", Ea=" << actual.lowPressure.Ea << std::endl;
                }

                // Verify Troe parameters
                if (!isEqual(actual.troe.a, 0.782, 0.01) ||
                    !isEqual(actual.troe.T_star, 207.5, 2.0) ||
                    !isEqual(actual.troe.T_double_star, 2663.0, 26.0) ||
                    !isEqual(actual.troe.T_triple_star, 6095.0, 60.0)) {
                    passed = false;
                    results.failureMessages.push_back(rxnEquation + ": Troe parameters mismatch");
                    std::cout << " " << rxnEquation << ": Troe parameters mismatch" << std::endl;
                    std::cout << "   Expected: a=0.782, T3=207.5, T1=2663.0, T2=6095.0" << std::endl;
                    std::cout << "   Actual: a=" << actual.troe.a << ", T3=" << actual.troe.T_star
                        << ", T1=" << actual.troe.T_double_star << ", T2=" << actual.troe.T_triple_star << std::endl;
                }

                if (passed) {
                    results.passedTests++;
                    std::cout << " Reaction " << rxnEquation << " test passed" << std::endl;
                }
            }
        }

        // Test results summary
        std::cout << "\n===== Specific C2H4 Data Test Summary =====\n";
        std::cout << "Total tests: " << results.totalTests << std::endl;
        std::cout << "Passed tests: " << results.passedTests << std::endl;
        std::cout << "Success rate: " << std::fixed << std::setprecision(2)
            << (results.totalTests > 0 ? 100.0 * results.passedTests / results.totalTests : 0)
            << "%" << std::endl;

        if (results.failureMessages.size() > 0) {
            std::cout << "\nFailed tests:\n";
            for (const auto& msg : results.failureMessages) {
                std::cout << " - " << msg << std::endl;
            }
        }

        return results.passedTests == results.totalTests;
    }
    catch (const std::exception& e) {
        std::cerr << "Error during specific C2H4 data test: " << e.what() << std::endl;
        return false;
    }
}


bool chemkinIOTest(const std::string& yamlFile)
{
    try {
        // Run specific C2H4 data test
        std::cout << "\n==============================================" << std::endl;
        std::cout << "      Specific C2H4 Data Test                " << std::endl;
        std::cout << "==============================================" << std::endl;
        bool specificTestResult = testSpecificC2H4Data(yamlFile);

        // Run general mechanism test
        std::cout << "==============================================" << std::endl;
        std::cout << "      General CH4 Data Test                  " << std::endl;
        std::cout << "==============================================" << std::endl;
        bool generalTestResult = runMechanismTest(yamlFile);

        // Combine results from all tests
        bool overallResult = generalTestResult && specificTestResult;

        // Print overall results
        if (overallResult) {
            std::cout << "\nOverall Test Result: All tests passed!" << std::endl;
        }
        else {
            std::cout << "\nOverall Test Result: Some tests failed. See details above." << std::endl;
        }

        return overallResult;
    }
    catch (const std::exception& e) {
        std::cerr << "Error during test execution: " << e.what() << std::endl;
        return false;
    }
}
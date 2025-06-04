# YAML Convector 2.0 - CMake Project

A chemical equilibrium calculation tool that processes YAML mechanism files and performs thermodynamic equilibrium calculations for ideal gas phases.

## Features

- Parse YAML chemistry mechanism files (compatible with Cantera format)
- Calculate chemical equilibrium for ideal gas mixtures
- Support for multiple thermodynamic property pairs (TP, HP, SP, etc.)
- Advanced Gibbs free energy minimization algorithms
- Cross-platform support (Windows, Linux, macOS)

## Requirements

### Build Requirements
- CMake 3.15 or higher
- C++14 compatible compiler
  - Visual Studio 2019/2022 (Windows)
  - GCC 7+ (Linux)
  - Clang 6+ (macOS)

### Dependencies
- yaml-cpp library (included in the project under `fluid.yaml-cpp.1.0.10/`)

## Quick Start

### Windows (Recommended)

Use the provided batch script for easy building:

```batch
# Configure the project (Release mode by default)
build.bat configure

# Build the project
build.bat build

# Run the executable
build.bat run

# Clean build files
build.bat clean
```

For Debug builds:
```batch
build.bat configure Debug
build.bat build
```

### Linux/macOS

Use the provided shell script:

```bash
# Make the script executable
chmod +x build.sh

# Configure the project
./build.sh configure

# Build the project
./build.sh build

# Run the executable
./build.sh run
```

### Manual CMake Commands (All Platforms)

```bash
# Create build directory
mkdir build
cd build

# Configure
cmake -DCMAKE_BUILD_TYPE=Release ..

# Build
cmake --build . --config Release

# Run (Windows)
./yaml-convector-2.0/Release/yaml-convector-2.0.exe

# Run (Linux/macOS)
./yaml-convector-2.0/yaml-convector-2.0
```

## Project Structure

```
yaml-convector-2.0/
├── CMakeLists.txt              # Main CMake configuration
├── build.bat                   # Windows build script
├── build.sh                    # Linux/macOS build script
├── README.md                   # This file
├── fluid.yaml-cpp.1.0.10/     # yaml-cpp library
│   └── build/
│       ├── include/            # Headers
│       └── lib/                # Libraries
├── yaml-convector-2.0/         # Source code directory
│   ├── CMakeLists.txt          # Source CMake configuration
│   ├── ChemEquil.cpp/.h        # Chemical equilibrium solver
│   ├── IdealGasPhase.cpp/.h    # Ideal gas phase implementation
│   ├── ChemistryIO.cpp/.h      # YAML file I/O
│   ├── ChemistryVars.cpp/.h    # Chemistry variables
│   ├── MechanismTest.cpp/.h    # Testing utilities
│   ├── yaml-convector-2.0.cpp # Main program
│   ├── ct_defs.h               # Type definitions
│   └── utils.h                 # Utility functions
└── build/                      # Build output (created by CMake)
```

## Usage Example

The program performs chemical equilibrium calculations for the H2/O2/Ar system:

```cpp
// Example from main program
IdealGasPhase gas;
gas.initFromYaml("D:\\cantera\\data\\h2o2.yaml", "ohmech");

// Set initial conditions: T=1500K, P=2atm, composition O2:1, H2:3, AR:1
gas.setState_TPX(1500.0, 2.0 * OneAtm, "O2:1.0, H2:3.0, AR:1.0");

// Perform equilibrium calculation
int result = gas.equilibrate("TP");
```

Expected output includes equilibrium compositions of all species (H2, H, O, O2, OH, H2O, HO2, H2O2, AR, N2).

## Build Configurations

### Debug Build
- Includes debug symbols
- No optimization
- Verbose error checking
- Links against yaml-cppd (debug version)

### Release Build (Default)
- Optimized for performance
- Minimal debug information
- Links against yaml-cpp (release version)

## Troubleshooting

### Common Issues

1. **CMake not found**
   - Install CMake from https://cmake.org/download/
   - Add CMake to your system PATH

2. **yaml-cpp library not found**
   - Ensure `fluid.yaml-cpp.1.0.10/build/` directory exists
   - Check that library files exist in `build/lib/`

3. **Build fails on Windows**
   - Use Visual Studio Developer Command Prompt
   - Ensure Visual Studio Build Tools are installed

4. **Executable not found after build**
   - Check build output directory
   - On Windows: `build/yaml-convector-2.0/Release/` or `build/yaml-convector-2.0/Debug/`
   - On Linux/macOS: `build/yaml-convector-2.0/`

### Debugging

To enable verbose CMake output:
```bash
cmake -DCMAKE_VERBOSE_MAKEFILE=ON ..
```

To see detailed build information:
```bash
cmake --build . --verbose
```

## Development

### Adding New Features
1. Add source files to `yaml-convector-2.0/CMakeLists.txt`
2. Update `YAML_CONVECTOR_SOURCES` and `YAML_CONVECTOR_HEADERS` lists
3. Rebuild the project

### Testing
The project includes built-in test functionality. Run with:
```bash
build.bat run    # Windows
./build.sh run   # Linux/macOS
```

## License

This project follows the same license terms as the original Cantera project. See individual source files for specific license information.

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Test thoroughly
5. Submit a pull request

For questions or issues, please create an issue in the project repository.

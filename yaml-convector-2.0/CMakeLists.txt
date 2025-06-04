cmake_minimum_required(VERSION 3.15)

# Project definition
project(yaml-convector-2.0 
    VERSION 2.0.0
    DESCRIPTION "YAML Convector for Chemical Equilibrium Calculations"
    LANGUAGES CXX
)

# Set C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set default build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Platform-specific settings
if(WIN32)
    add_definitions(-DWIN32 -D_CONSOLE)
    if(MSVC)
        add_compile_options(/W3)
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MDd")
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MD")
    endif()
else()
    add_compile_options(-Wall -Wextra)
endif()

# Find required packages
find_package(PkgConfig)

# Try to find yaml-cpp
find_package(yaml-cpp QUIET)

if(NOT yaml-cpp_FOUND)
    # If yaml-cpp is not found via find_package, try to find it manually
    set(YAML_CPP_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/fluid.yaml-cpp.1.0.10")
    
    if(EXISTS "${YAML_CPP_ROOT}/build/include/yaml-cpp/yaml.h")
        message(STATUS "Found local yaml-cpp installation at ${YAML_CPP_ROOT}")
        
        # Create imported target for yaml-cpp
        add_library(yaml-cpp SHARED IMPORTED)
        
        # Set include directories
        set_target_properties(yaml-cpp PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${YAML_CPP_ROOT}/build/include"
        )
        
        # Set library location based on build type and platform
        if(WIN32)
            if(CMAKE_BUILD_TYPE STREQUAL "Debug")
                set_target_properties(yaml-cpp PROPERTIES
                    IMPORTED_LOCATION "${YAML_CPP_ROOT}/build/lib/yaml-cppd.dll"
                    IMPORTED_IMPLIB "${YAML_CPP_ROOT}/build/lib/yaml-cppd.lib"
                )
            else()
                set_target_properties(yaml-cpp PROPERTIES
                    IMPORTED_LOCATION "${YAML_CPP_ROOT}/build/lib/yaml-cpp.dll"
                    IMPORTED_IMPLIB "${YAML_CPP_ROOT}/build/lib/yaml-cpp.lib"
                )
            endif()
        else()
            # For non-Windows platforms, assume static library
            if(CMAKE_BUILD_TYPE STREQUAL "Debug")
                set_target_properties(yaml-cpp PROPERTIES
                    IMPORTED_LOCATION "${YAML_CPP_ROOT}/build/lib/libyaml-cppd.a"
                )
            else()
                set_target_properties(yaml-cpp PROPERTIES
                    IMPORTED_LOCATION "${YAML_CPP_ROOT}/build/lib/libyaml-cpp.a"
                )
            endif()
        endif()
        
        set(yaml-cpp_FOUND TRUE)
    else()
        message(WARNING "yaml-cpp not found. Please install yaml-cpp or provide the path.")
    endif()
endif()

# Add the main executable target
add_subdirectory(yaml-convector-2.0)

# Configure installation
include(GNUInstallDirs)

# Install the executable
install(TARGETS yaml-convector-2.0-exe 
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# Install yaml files
install(FILES 
    "${CMAKE_CURRENT_SOURCE_DIR}/C2H4-mechanism.yaml"
    "${CMAKE_CURRENT_SOURCE_DIR}/mechanism.yaml"
    DESTINATION ${CMAKE_INSTALL_DATADIR}/yaml-convector-2.0
    OPTIONAL
)

# Copy DLLs on Windows for runtime
if(WIN32 AND yaml-cpp_FOUND)
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        install(FILES "${YAML_CPP_ROOT}/build/lib/yaml-cppd.dll"
            DESTINATION ${CMAKE_INSTALL_BINDIR}
            OPTIONAL
        )
    else()
        install(FILES "${YAML_CPP_ROOT}/build/lib/yaml-cpp.dll"
            DESTINATION ${CMAKE_INSTALL_BINDIR}
            OPTIONAL
        )
    endif()
endif()

# Print configuration summary
message(STATUS "")
message(STATUS "=== Configuration Summary ===")
message(STATUS "Project: ${PROJECT_NAME} v${PROJECT_VERSION}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "yaml-cpp found: ${yaml-cpp_FOUND}")
if(yaml-cpp_FOUND AND TARGET yaml-cpp)
    get_target_property(YAML_CPP_INCLUDE yaml-cpp INTERFACE_INCLUDE_DIRECTORIES)
    message(STATUS "yaml-cpp include: ${YAML_CPP_INCLUDE}")
endif()
message(STATUS "=============================")
message(STATUS "")

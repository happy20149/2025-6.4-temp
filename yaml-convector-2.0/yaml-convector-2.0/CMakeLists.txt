# Source directory CMakeLists.txt

# Define source files
set(YAML_CONVECTOR_SOURCES
    ChemEquil.cpp
    ChemistryIO.cpp
    ChemistryVars.cpp
    IdealGasPhase.cpp
    MechanismTest.cpp
    yaml-convector-2.0.cpp
)

# Define header files
set(YAML_CONVECTOR_HEADERS
    ChemEquil.h
    ChemistryIO.h
    ChemistryVars.h
    IdealGasPhase.h
    MechanismTest.h
    ct_defs.h
    utils.h
)

# Create the executable
add_executable(yaml-convector-2.0-exe ${YAML_CONVECTOR_SOURCES})

# Set target properties
set_target_properties(yaml-convector-2.0-exe PROPERTIES
    OUTPUT_NAME yaml-convector-2.0
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
)

# Include directories
target_include_directories(yaml-convector-2.0-exe 
    PRIVATE 
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Link libraries
if(TARGET yaml-cpp)
    target_link_libraries(yaml-convector-2.0-exe PRIVATE yaml-cpp)
    
    # Add yaml-cpp definitions if needed
    target_compile_definitions(yaml-convector-2.0-exe PRIVATE 
        HAVE_YAML_CPP     # Enable yaml-cpp support
        YAML_CPP_DLL      # For Windows DLL
    )
else()
    message(WARNING "yaml-cpp not found. Building without yaml-cpp support.")
    target_compile_definitions(yaml-convector-2.0-exe PRIVATE 
        NO_YAML_CPP      # Disable yaml-cpp support
    )
endif()

# Platform-specific settings
if(WIN32)
    # Windows-specific settings
    target_compile_definitions(yaml-convector-2.0-exe PRIVATE
        WIN32
        _CONSOLE
        UNICODE
        _UNICODE
    )
    
    # Set subsystem to console
    set_target_properties(yaml-convector-2.0-exe PROPERTIES
        WIN32_EXECUTABLE FALSE
    )
    
    if(MSVC)
        # MSVC-specific settings
        target_compile_options(yaml-convector-2.0-exe PRIVATE
            /W3          # Warning level 3
            /permissive- # Conformance mode
        )
        
        # Runtime library settings
        if(CMAKE_BUILD_TYPE STREQUAL "Debug")
            target_compile_options(yaml-convector-2.0-exe PRIVATE /MDd)
            target_compile_definitions(yaml-convector-2.0-exe PRIVATE _DEBUG)
        else()
            target_compile_options(yaml-convector-2.0-exe PRIVATE /MD)
            target_compile_definitions(yaml-convector-2.0-exe PRIVATE NDEBUG)
        endif()
    endif()
else()
    # Unix-like systems settings
    target_compile_options(yaml-convector-2.0-exe PRIVATE
        -Wall
        -Wextra
        -Wpedantic
    )
    
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        target_compile_options(yaml-convector-2.0-exe PRIVATE -g -O0)
        target_compile_definitions(yaml-convector-2.0-exe PRIVATE DEBUG)
    else()
        target_compile_options(yaml-convector-2.0-exe PRIVATE -O3)
        target_compile_definitions(yaml-convector-2.0-exe PRIVATE NDEBUG)
    endif()
endif()

# Add custom target for running the executable
add_custom_target(run
    COMMAND yaml-convector-2.0-exe
    DEPENDS yaml-convector-2.0-exe
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Running yaml-convector-2.0"
)

# Add test target if needed
# add_custom_target(test
#     COMMAND yaml-convector-2.0-exe --test
#     DEPENDS yaml-convector-2.0-exe
#     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
#     COMMENT "Running yaml-convector-2.0 tests"
# )

# Print target information
message(STATUS "Adding executable: yaml-convector-2.0-exe")
message(STATUS "Source files: ${YAML_CONVECTOR_SOURCES}")
message(STATUS "Header files: ${YAML_CONVECTOR_HEADERS}")

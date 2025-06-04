#!/bin/bash

# Build script for yaml-convector-2.0 CMake project
# This script provides easy commands to configure, build, and run the project

set -e

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="$PROJECT_ROOT/build"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

print_usage() {
    echo "Usage: $0 [command] [options]"
    echo ""
    echo "Commands:"
    echo "  configure [Release|Debug]  - Configure the build (default: Release)"
    echo "  build                      - Build the project"
    echo "  clean                      - Clean build directory"
    echo "  rebuild                    - Clean and build"
    echo "  run                        - Run the executable"
    echo "  install                    - Install the project"
    echo "  help                       - Show this help"
    echo ""
    echo "Examples:"
    echo "  $0 configure Debug"
    echo "  $0 build"
    echo "  $0 run"
}

print_status() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

configure_project() {
    local build_type=${1:-Release}
    
    print_status "Configuring project with build type: $build_type"
    
    # Create build directory
    mkdir -p "$BUILD_DIR"
    cd "$BUILD_DIR"
    
    # Configure with CMake
    cmake \
        -DCMAKE_BUILD_TYPE="$build_type" \
        -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
        "$PROJECT_ROOT"
    
    print_success "Configuration completed successfully"
}

build_project() {
    if [ ! -d "$BUILD_DIR" ]; then
        print_error "Build directory not found. Please run 'configure' first."
        exit 1
    fi
    
    print_status "Building project..."
    cd "$BUILD_DIR"
    
    # Build with CMake
    cmake --build . --config Release
    
    print_success "Build completed successfully"
    print_status "Executable location: $BUILD_DIR/yaml-convector-2.0/yaml-convector-2.0"
}

clean_project() {
    print_status "Cleaning build directory..."
    if [ -d "$BUILD_DIR" ]; then
        rm -rf "$BUILD_DIR"
        print_success "Clean completed"
    else
        print_warning "Build directory doesn't exist"
    fi
}

run_project() {
    if [ ! -d "$BUILD_DIR" ]; then
        print_error "Build directory not found. Please build the project first."
        exit 1
    fi
    
    local executable="$BUILD_DIR/yaml-convector-2.0/yaml-convector-2.0"
    if [ ! -f "$executable" ]; then
        print_error "Executable not found. Please build the project first."
        exit 1
    fi
    
    print_status "Running yaml-convector-2.0..."
    cd "$BUILD_DIR/yaml-convector-2.0"
    ./yaml-convector-2.0
}

install_project() {
    if [ ! -d "$BUILD_DIR" ]; then
        print_error "Build directory not found. Please build the project first."
        exit 1
    fi
    
    print_status "Installing project..."
    cd "$BUILD_DIR"
    cmake --install .
    print_success "Installation completed"
}

# Main script logic
case "${1:-help}" in
    configure)
        configure_project "$2"
        ;;
    build)
        build_project
        ;;
    clean)
        clean_project
        ;;
    rebuild)
        clean_project
        configure_project "$2"
        build_project
        ;;
    run)
        run_project
        ;;
    install)
        install_project
        ;;
    help)
        print_usage
        ;;
    *)
        print_error "Unknown command: $1"
        print_usage
        exit 1
        ;;
esac

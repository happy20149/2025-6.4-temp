@echo off
setlocal enabledelayedexpansion

:: Build script for yaml-convector-2.0 CMake project (Windows)
:: This script provides easy commands to configure, build, and run the project

set "PROJECT_ROOT=%~dp0"
set "BUILD_DIR=%PROJECT_ROOT%build"

goto :main

:print_usage
echo Usage: %~nx0 [command] [options]
echo.
echo Commands:
echo   configure [Release^|Debug]  - Configure the build (default: Release)
echo   build                      - Build the project
echo   clean                      - Clean build directory
echo   rebuild                    - Clean and build
echo   run                        - Run the executable
echo   install                    - Install the project
echo   help                       - Show this help
echo.
echo Examples:
echo   %~nx0 configure Debug
echo   %~nx0 build
echo   %~nx0 run
goto :eof

:print_status
echo [INFO] %~1
goto :eof

:print_success
echo [SUCCESS] %~1
goto :eof

:print_warning
echo [WARNING] %~1
goto :eof

:print_error
echo [ERROR] %~1
goto :eof

:configure_project
set "build_type=%~1"
if "%build_type%"=="" set "build_type=Release"

call :print_status "Configuring project with build type: %build_type%"

:: Create build directory
if not exist "%BUILD_DIR%" mkdir "%BUILD_DIR%"
cd /d "%BUILD_DIR%"

:: Try to find CMake
where cmake >nul 2>&1
if errorlevel 1 (
    call :print_error "CMake not found in PATH. Please install CMake or add it to PATH."
    exit /b 1
)

:: Configure with CMake
cmake ^
    -DCMAKE_BUILD_TYPE="%build_type%" ^
    -DCMAKE_EXPORT_COMPILE_COMMANDS=ON ^
    "%PROJECT_ROOT%"

if errorlevel 1 (
    call :print_error "Configuration failed"
    exit /b 1
)

call :print_success "Configuration completed successfully"
goto :eof

:build_project
if not exist "%BUILD_DIR%" (
    call :print_error "Build directory not found. Please run 'configure' first."
    exit /b 1
)

call :print_status "Building project..."
cd /d "%BUILD_DIR%"

:: Build with CMake
cmake --build . --config Release

if errorlevel 1 (
    call :print_error "Build failed"
    exit /b 1
)

call :print_success "Build completed successfully"
call :print_status "Executable location: %BUILD_DIR%\yaml-convector-2.0\Release\yaml-convector-2.0.exe"
goto :eof

:clean_project
call :print_status "Cleaning build directory..."
if exist "%BUILD_DIR%" (
    rmdir /s /q "%BUILD_DIR%"
    call :print_success "Clean completed"
) else (
    call :print_warning "Build directory doesn't exist"
)
goto :eof

:run_project
if not exist "%BUILD_DIR%" (
    call :print_error "Build directory not found. Please build the project first."
    exit /b 1
)

:: Try different possible executable locations
set "executable1=%BUILD_DIR%\yaml-convector-2.0\Release\yaml-convector-2.0.exe"
set "executable2=%BUILD_DIR%\yaml-convector-2.0\Debug\yaml-convector-2.0.exe"
set "executable3=%BUILD_DIR%\yaml-convector-2.0\yaml-convector-2.0.exe"

set "executable="
if exist "%executable1%" set "executable=%executable1%"
if exist "%executable2%" set "executable=%executable2%"
if exist "%executable3%" set "executable=%executable3%"

if "%executable%"=="" (
    call :print_error "Executable not found. Please build the project first."
    call :print_status "Looked for:"
    echo   %executable1%
    echo   %executable2%
    echo   %executable3%
    exit /b 1
)

call :print_status "Running yaml-convector-2.0..."
cd /d "%~dp1"
"%executable%"
goto :eof

:install_project
if not exist "%BUILD_DIR%" (
    call :print_error "Build directory not found. Please build the project first."
    exit /b 1
)

call :print_status "Installing project..."
cd /d "%BUILD_DIR%"
cmake --install .

if errorlevel 1 (
    call :print_error "Installation failed"
    exit /b 1
)

call :print_success "Installation completed"
goto :eof

:main
set "command=%~1"
if "%command%"=="" set "command=help"

if "%command%"=="configure" (
    call :configure_project "%~2"
) else if "%command%"=="build" (
    call :build_project
) else if "%command%"=="clean" (
    call :clean_project
) else if "%command%"=="rebuild" (
    call :clean_project
    call :configure_project "%~2"
    call :build_project
) else if "%command%"=="run" (
    call :run_project
) else if "%command%"=="install" (
    call :install_project
) else if "%command%"=="help" (
    call :print_usage
) else (
    call :print_error "Unknown command: %command%"
    call :print_usage
    exit /b 1
)

exit /b 0

@echo off
echo Configuring CMake project...
if not exist "build" mkdir build
cd build
cmake .. -G "Visual Studio 17 2022" -A x64 -DCMAKE_TOOLCHAIN_FILE=C:/Users/Administrator/vcpkg/scripts/buildsystems/vcpkg.cmake
if %errorlevel% neq 0 (
    echo CMake configuration failed!
    pause
    exit /b %errorlevel%
)

echo Building project...
cmake --build . --config Release
if %errorlevel% neq 0 (
    echo Build failed!
    pause
    exit /b %errorlevel%
)

echo Build completed successfully!
pause

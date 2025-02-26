:: bld.bat
:: For Windows builds (Not tested)
@echo off
setlocal EnableDelayedExpansion

mkdir build
cd build

:: Check if Ninja is available
where ninja >nul 2>nul
if %ERRORLEVEL% EQU 0 (
    echo Using Ninja
    set GENERATOR=Ninja
) else (
    echo Using NMake Makefiles
    set GENERATOR=NMake Makefiles
)

:: Run CMake Configuration
cmake -G "!GENERATOR!" ^
      -DCMAKE_CXX_COMPILER="%CXX%" ^
      -DCMAKE_C_FLAGS="/DHAVE_STRDUP=1" ^
      -DCMAKE_INSTALL_PREFIX="%PREFIX%" ^
      ..
if %ERRORLEVEL% NEQ 0 (
    echo ERROR: CMake configuration failed! >&2
    exit /b %ERRORLEVEL%
)

:: Build
cmake --build . --parallel
if %ERRORLEVEL% NEQ 0 (
    echo ERROR: Build failed! >&2
    exit /b %ERRORLEVEL%
)

:: Install
cmake --install .
if %ERRORLEVEL% NEQ 0 (
    echo ERROR: Installation failed! >&2
    exit /b %ERRORLEVEL%
)

echo Build and installation completed successfully.

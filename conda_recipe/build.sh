#!/bin/bash
mkdir -p build
cd build

# Check if Ninja is installed
if command -v ninja &>/dev/null; then
    echo "Using Ninja"
    GENERATOR="Ninja"
else
    echo "Using Unix Makefiles"
    GENERATOR="Unix Makefiles"
fi

# Run CMake with error handling
cmake -G "$GENERATOR" \
      -DCMAKE_CXX_COMPILER="${CXX:-g++}" \
      -DCMAKE_C_FLAGS="-DHAVE_STRDUP=1" \
      -DCMAKE_INSTALL_PREFIX="${PREFIX:-/usr/local}" \
      .. || { echo "CMake configuration failed"; exit 1; }

# Build
cmake --build . --parallel || { echo "Build failed"; exit 1; }

# Install
cmake --install . || { echo "Installation failed"; exit 1; }

echo "Build and installation completed successfully."

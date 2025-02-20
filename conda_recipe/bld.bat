:: Not tested (For windows build)
mkdir build
cd build
cmake -G "NMake Makefiles" -DCMAKE_CXX_COMPILER=%CXX% -DCMAKE_C_FLAGS="-DHAVE_STRDUP=1" ..
cmake --build . --parallel $(nproc)
cmake --install . --prefix=%PREFIX%

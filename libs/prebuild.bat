@echo "Building libdeflate"
cd libdeflate

cmake -B build
cmake --build build --config Debug
cmake --build build --config Release


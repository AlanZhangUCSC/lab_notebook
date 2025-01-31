startDir=$PWD  # Note: should be PWD not pwd
cd $(dirname "$0")
mkdir -p ../build
cd ../build

# Install capnp with proper prefix and pkg-config setup
curl -O https://capnproto.org/capnproto-c++-1.0.2.tar.gz
tar zxf capnproto-c++-1.0.2.tar.gz
cd capnproto-c++-1.0.2
./configure --prefix=/usr
make -j check
make install
ldconfig  # Update shared library cache
cd ../

# Clone and setup vcpkg
git clone https://github.com/microsoft/vcpkg.git
apt-get install pkg-config
./vcpkg/bootstrap-vcpkg.sh
./vcpkg/vcpkg install jsoncpp

# Download and extract TBB
wget https://github.com/oneapi-src/oneTBB/archive/2019_U9.tar.gz
tar -xvzf 2019_U9.tar.gz

# Run CMake with proper paths
cmake \
    -DTBB_DIR=${PWD}/oneTBB-2019_U9 \
    -DCMAKE_PREFIX_PATH="${PWD}/oneTBB-2019_U9/cmake;/usr" \
    -DCMAKE_TOOLCHAIN_FILE=${PWD}/vcpkg/scripts/buildsystems/vcpkg.cmake \
    ..

make -j
cd $startDir
startDir=$PWD  # Note: should be PWD not pwd
cd $(dirname "$0")
mkdir -p ../build
cd ../build


cd capnproto-c++-1.0.2
./configure --prefix=/usr
make -j check
make install
ldconfig  # Update shared library cache
cd ../

cd $startDir
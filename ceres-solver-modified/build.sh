sudo rm -r build
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=./ .. -DCERES_THREADING_MODEL=OPENMP
make -j8
sudo make install


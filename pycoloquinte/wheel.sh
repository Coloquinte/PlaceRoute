# Run on docker run -it quay.io/pypa/manylinux2014_x86_64 /bin/bash

# Install the dependencies
yum install boost-devel boost-static eigen3-devel -y
# Download and install lemon; package coin-or-lemon-devel not available
yum install wget -y
wget http://lemon.cs.elte.hu/pub/sources/lemon-1.3.1.tar.gz; tar -xzf lemon-1.3.1.tar.gz; cd lemon-1.3.1; mkdir build; cd build; cmake ..; make install; cd ../..

# Build the wheels
git clone https://github.com/Coloquinte/PlaceRoute --recursive
cd PlaceRoute/pycoloquinte/
for p in /usr/local/bin/python3*; do ${p} setup.py bdist_wheel; done
for w in dist/coloquinte*; do auditwheel repair $w; done

# Upload the wheels
python3.11 -m pip install twine
python3.11 -m twine upload wheelhouse/*

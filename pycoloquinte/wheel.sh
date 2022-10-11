# Run on docker run -it quay.io/pypa/manylinux2014_x86_64 /bin/bash
git clone https://github.com/Coloquinte/PlaceRoute --recursive
yum install boost-devel boost-static eigen3-devel -y
cd PlaceRoute/pycoloquinte/
for p in /usr/local/bin/python3*; do ${p} setup.py bdist_wheel; done
for w in dist/coloquinte*; do auditwheel repair $w; done
python3.11 -m pip install twine
python3.11 -m twine upload wheelhouse/*

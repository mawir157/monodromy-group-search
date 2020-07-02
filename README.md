
# monodromy-group-search

Program for exploring the properties of complex hyperbolic triangle groups this is an improved and
cannibalized version of complex-hyperbolic-triangle-explorer.

## Installation
Installation requires the Armadillo linear algebra library and cmake. 

1. Install OpenSSL `sudo apt-get install libssl-dev`

2. Install Blas and LAAPACK `sudo apt-get install libblas-dev liblapack-dev`

2. Install cmake https://cmake.org/install/

3. Install Armadillo

wget http://sourceforge.net/projects/arma/files/armadillo-x.yyy.z.tar.gz
or go to http://arma.sourceforge.net/download.html
```
tar -xvf armadillo-x.yyy.z.tar.gz
cd armadillo-x.yyy.z
./configure
make
sudo make install
```

5. cd into monodromy-group-search then 'make'

## Usage

`mono -f input/groups.txt`

`mono -l 5 5`

`mono -m input/matrices.txt`

### Mode -m
This mode allows the user to input their own matrices in a file (see ./input/matrices.txt for an example). The code uses a tolerance of between 1e-6 and 1e-10 to check when two complex numbers are equal. Consequently the values in the text file should be given to around 15 decimal places otherwise rounding errors start to creep in and code fails in unpredictable ways where is difficult to tell if something is real or an artefact of rounding.

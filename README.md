# PoKiTT 
Portable Kinetics Thermodynamics & Transport (PoKiTT) is a library that
facilitates calculation of a variety of quantities related to reacting and
nonreacting flows.

------------------------------

# Building PoKiTT 
## Downloading the source code
The source may be obtained via git:
```
git clone https://gitlab.multiscale.utah.edu/common/PoKiTT.git [destination path]
```


## Cantera dependency

PoKiTT depends on [Cantera](www.cantera.org) to parse input files and generate data structures.
PoKiTT interrogates these data structures to build its own, optimized, calculations of thermodymic and transport properties as well as reaction rates and sensitivities.

We have a [modified version of Cantera](https://gitlab.multiscale.utah.edu/common/cantera.git) which provides additional API functionality that is required for PoKiTT.
Before installing PoKiTT, first install [that version of Cantera](https://gitlab.multiscale.utah.edu/common/cantera.git).


## Building PoKiTT for Downstream apps

```
git clone https://gitlab.multiscale.utah.edu/common/PoKiTT.git
cd PoKiTT; mkdir build; cd build
cmake .. \
  -DBUILD_UPSTREAM_LIBS=ON \
  -DCantera_DIR=[path to cantera install] \
  -DExprLib_DIR=[path to exprlib install] \
  -DENABLE_TESTS=OFF
```


## Building PoKiTT for Stand-Alone use (testing, etc)

PoKiTT relies on two upstream libraries: 
  * [`SpatialOps`](https://gitlab.multiscale.utah.edu/common/SpatialOps) provides low-level data parallelism over grid operations and supports both CPU and GPU execution modes.
  * [`ExprLib`](https://gitlab.multiscale.utah.edu/common/ExprLib) provides task parallelism and manages complex data and execution dependencies in multiphysics software.

### A simple CPU-Only Build

The simplest way to build the examples and tests is to automatically build [`SpatialOps`](https://gitlab.multiscale.utah.edu/common/SpatialOps) and [`ExprLib`](https://gitlab.multiscale.utah.edu/common/ExprLib) by setting the `BUILD_UPSTREAM_LIBS` flag to `ON`:
```
git clone https://gitlab.multiscale.utah.edu/common/PoKiTT.git
cd PoKiTT; mkdir build; cd build
cmake .. \
  -DBUILD_UPSTREAM_LIBS=ON
  -DCantera_DIR=[path to cantera install] \
  -DENABLE_TESTS=ON
```
This produces a serial build of PoKiTT; no multithreaded or GPU capabilities are active.


### A GPU Build

This assumes that you are building in a directory defined by the variable
`BUILD_DIR` and that the source is in a directory defined by the variable `SRC_DIR`
```
cd $SRC_DIR
git clone https://gitlab.multiscale.utah.edu/common/PoKiTT.git .

cd $BUILD_DIR

mkdir tpl; cd tpl;
git clone --depth=1 https://gitlab.multiscale.utah.edu/common/SpatialOps.git src/SpatialOps
git clone --depth=1 https://gitlab.multiscale.utah.edu/common/ExprLib.git src/ExprLib

mkdir so; cd so;
cmake $BUILD_DIR/tpl/src/SpatialOps -DENABLE_CUDA=ON -DENABLE_TESTS=OFF -DENABLE_THREADS=OFF -DENABLE_EXAMPLES=OFF -DCMAKE_INSTALL_PREFIX=${BUILD_DIR}/tpl
make -j8 install

cd $BUILD_DIR/tpl
mkdir expr; cd expr;
cmake $BUILD_DIR/tpl/src/ExprLib -DENABLE_TESTS=OFF -DCMAKE_INSTALL_PREFIX=$BUILD_DIR/tpl -DSpatialOps_DIR=$BUILD_DIR/tpl/share
make -j8 install

cd $BUILD_DIR
cmake $SRC_DIR -DExprLib_DIR=$BUILD_DIR/tpl/share -DCantera_DIR=[path to Cantera]
make -j8 install
```

#!/bin/bash

if [[ $# != 1 ]]; then
    echo "Correct format: compile.sh <ALGORITHM>"
    exit 1
fi

WORK_DIR="$(pwd)"
cd "$WORK_DIR/bin/.o/"

if [[ $1 == "gorilla" ]]; then
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++20 -c "$WORK_DIR"/src/c++/floating-point/lossless/gorilla.cpp
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++20 -c "$WORK_DIR"/src/c++/main.cpp
elif [[ $1 == "chimp" ]]; then
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++20 -c "$WORK_DIR"/src/c++/floating-point/lossless/chimp.cpp
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++20 -c "$WORK_DIR"/src/c++/main.cpp
elif [[ $1 == "elf" ]]; then
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++20 -c "$WORK_DIR"/src/c++/floating-point/lossless/elf.cpp
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++20 -c "$WORK_DIR"/src/c++/main.cpp
elif [[ $1 == "self" ]]; then
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++20 -c "$WORK_DIR"/src/c++/floating-point/lossless/self.cpp
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++20 -c "$WORK_DIR"/src/c++/main.cpp
elif [[ $1 == "serf" ]]; then
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++20 -c "$WORK_DIR"/src/c++/floating-point/lossy/serf.cpp
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++20 -c "$WORK_DIR"/src/c++/main.cpp
elif [[ $1 == "camel" ]]; then
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++20 -c "$WORK_DIR"/src/c++/floating-point/lossy/camel.cpp
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++20 -c "$WORK_DIR"/src/c++/main.cpp
elif [[ $1 == "dexor" ]]; then
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++20  -c "$WORK_DIR"/src/c++/floating-point/lossless/dexor.cpp
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++20 -c "$WORK_DIR"/src/c++/main.cpp
elif [[ $1 == "kangaroo" ]]; then
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++20  -c "$WORK_DIR"/src/c++/floating-point/lossless/kangaroo.cpp
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++20 -c "$WORK_DIR"/src/c++/main.cpp
elif [[ $1 == "pmc" ]]; then
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++20 -c "$WORK_DIR"/src/c++/piecewise-approximation/constant/pmc.cpp
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++20 -c "$WORK_DIR"/src/c++/main.cpp
elif [[ $1 == "hybrid-pca" ]]; then
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++20 -c "$WORK_DIR"/src/c++/piecewise-approximation/constant/hybrid-pca.cpp
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++20 -c "$WORK_DIR"/src/c++/main.cpp
elif [[ $1 == "swing-filter" ]]; then
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++20 -c "$WORK_DIR"/src/c++/piecewise-approximation/linear/swing-filter.cpp
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++20 -c "$WORK_DIR"/src/c++/main.cpp
elif [[ $1 == "slide-filter" ]]; then
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++20 -c "$WORK_DIR"/src/c++/piecewise-approximation/linear/slide-filter.cpp
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++20 -c "$WORK_DIR"/src/c++/main.cpp
elif [[ $1 == "optimal-pla" ]]; then
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++20 -c "$WORK_DIR"/src/c++/piecewise-approximation/linear/optimal-pla.cpp
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++20 -c "$WORK_DIR"/src/c++/main.cpp
elif [[ $1 == "cov-pla" ]]; then
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++20 -c "$WORK_DIR"/src/c++/piecewise-approximation/linear/cov-pla.cpp
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++20 -c "$WORK_DIR"/src/c++/main.cpp
elif [[ $1 == "conn-I-pla" ]]; then
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++20 -c "$WORK_DIR"/src/c++/piecewise-approximation/linear/conn-I-pla.cpp
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++20 -c "$WORK_DIR"/src/c++/main.cpp
elif [[ $1 == "semi-optimal-pla" ]]; then
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++20 -c "$WORK_DIR"/src/c++/piecewise-approximation/linear/semi-optimal-pla.cpp
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++20 -c "$WORK_DIR"/src/c++/main.cpp
elif [[ $1 == "semi-mixed-pla" ]]; then
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++20 -c "$WORK_DIR"/src/c++/piecewise-approximation/linear/semi-mixed-pla.cpp
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++20 -c "$WORK_DIR"/src/c++/main.cpp
elif [[ $1 == "mix-piece" ]]; then
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++20 -c "$WORK_DIR"/src/c++/piecewise-approximation/linear/mix-piece.cpp
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++20 -c "$WORK_DIR"/src/c++/main.cpp
elif [[ $1 == "ioriented-pla" ]]; then
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++20 -c "$WORK_DIR"/src/c++/piecewise-approximation/linear/ioriented-pla.cpp
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++20 -c "$WORK_DIR"/src/c++/main.cpp
elif [[ $1 == "cached-normal-equation" ]]; then
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++20 -c "$WORK_DIR"/src/c++/piecewise-approximation/polynomial/cached-normal-equation.cpp
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++20 -c "$WORK_DIR"/src/c++/main.cpp
elif [[ $1 == "poly-swab" ]]; then
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++20 -c "$WORK_DIR"/src/c++/piecewise-approximation/polynomial/poly-swab.cpp
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++20 -c "$WORK_DIR"/src/c++/main.cpp
elif [[ $1 == "swab" ]]; then
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++20 -c "$WORK_DIR"/src/c++/piecewise-approximation/linear/swab.cpp
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++20 -c "$WORK_DIR"/src/c++/main.cpp
elif [[ $1 == "addaptive-approximation" ]]; then
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++20 -c "$WORK_DIR"/src/c++/model-selection/adaptive-approximation.cpp
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++20 -c "$WORK_DIR"/src/c++/main.cpp
elif [[ $1 == "smart-grid-compression" ]]; then
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++20 -c "$WORK_DIR"/src/c++/model-selection/smart-grid-compression.cpp
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++20 -c "$WORK_DIR"/src/c++/main.cpp
elif [[ $1 == "adapt-ppa" ]]; then
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++20 -c "$WORK_DIR"/src/c++/model-selection/adapt-ppa.cpp
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++20 -c "$WORK_DIR"/src/c++/main.cpp
else
    echo "Algorithm $1 is not supported"
    exit 2
fi

cd "$WORK_DIR"/bin/
g++ $(find $WORK_DIR/bin/.o/ -name "*.o" ! -name "main.o" -type f | xargs) "$WORK_DIR"/bin/.o/main.o -o main

exit 0

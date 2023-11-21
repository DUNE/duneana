#!/bin/bash

echo "COMPILING CLUSTERING MODULE."
g++ -std=c++11 -o Module_SNClustering.exe Module_SNClustering.C `root-config --cflags --glibs`
echo "COMPILING CLUSTERING ANALYSIS."
g++ -std=c++11 -o Analyse_SNClustering.exe Analyse_SNClustering.C `root-config --cflags --glibs`

#!/bin/sh
set -e
# Build.sh
# mesmer
#
# Created by Struan Robertson on 13/03/2012.


echo "----------------------------------------------------------------------"
echo " Building TinyXML."
echo "----------------------------------------------------------------------"

cd ./tinyxml
make -j 4 -f MakeLib DEBUG=NO
cd ../

echo "----------------------------------------------------------------------"
echo " Building QD."
echo "----------------------------------------------------------------------"

cd ./qd
chmod +x configure
CXX=gcc CXXFLAGS='-O2' ./configure --prefix $PWD/src
make -j 4 install
cd ../

echo "----------------------------------------------------------------------"
echo " Building Mesmer."
echo "----------------------------------------------------------------------"

cd ./src
make -j 4 install DEBUG=NO
cd ../


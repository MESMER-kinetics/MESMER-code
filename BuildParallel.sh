#!/bin/sh

# Build.sh
# mesmer
#
# Created by Struan Robertson on 13/03/2012.


echo "----------------------------------------------------------------------"
echo " Building TinyXML."
echo "----------------------------------------------------------------------"

cd ./tinyxml
make -f MakeLib DEBUG=NO
cd ../

echo "----------------------------------------------------------------------"
echo " Building QD."
echo "----------------------------------------------------------------------"

cd ./qd
chmod +x configure
./configure
make
cd ../

echo "----------------------------------------------------------------------"
echo " Building Mesmer."
echo "----------------------------------------------------------------------"

cd ./src
make -f Makefile install DEBUG=NO PARALLEL=YES
cd ../


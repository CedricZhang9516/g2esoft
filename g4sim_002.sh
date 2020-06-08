#!/bin/bash

DIR=/pbs/home/t/tyamanak/private/sps/tyamanak/g2esoft/yamanaka/g2esoft

source /pbs/throng/g-2/share/setup_gcc9_2_0.sh
export LD_LIBRARY_PATH=$DIR/lib:${LD_LIBRARY_PATH}
cp $DIR/config/testG4Sim.xml .
sed -e 's/<int name="EventNumToProcess">100<\/int>/<int name="EventNumToProcess">500<\/int>/g' testG4Sim.xml > temp1.xml
sed -e 's/<int name="EventIdxToStart">0<\/int>/<int name="EventIdxToStart">500<\/int>/g' temp1.xml > temp2.xml
sed -e 's/<int name="RandomSeed">1<\/int>/<int name="RandomSeed">2<\/int>/g' temp2.xml > testG4Sim_final.xml
cmake $DIR
make install
./g2esoft.bin testG4Sim_final.xml
cp testG4Out.root ${HOME}/G4Out_002.root

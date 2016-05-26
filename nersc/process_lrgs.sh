#!/bin/bash
CONF=/project/projectdirs/m779/yfeng1/imaginglss/dr2.conf.py

cp RANDOM.hdf5 LRG-RANDOM.hdf5

python ../scripts/imglss-query-tycho-veto.py LRG-RANDOM.hdf5 --conf $CONF
python ../scripts/imglss-query-tycho-veto.py LRG.hdf5 --conf $CONF
python ../scripts/imglss-query-completeness.py --conf $CONF \
        --use-tycho-veto=BOSS_DR9 LRG LRG.hdf5 LRG.hdf5
python ../scripts/imglss-query-completeness.py --conf $CONF \
        --use-tycho-veto=BOSS_DR9 LRG LRG.hdf5 LRG-RANDOM.hdf5
python '../scripts/imglss-export-text.py' --conf $CONF --sigma-g -1000 \
        --use-tycho-veto=BOSS_DR9 --bands r g z W1 W2 -- LRG.hdf5 LRG.txt
python '../scripts/imglss-export-text.py' --conf $CONF --sigma-g -1000 \
        --use-tycho-veto=BOSS_DR9 --bands r g z W1 W2 -- LRG-RANDOM.hdf5 LRG-RANDOM.txt
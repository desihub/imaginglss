#!/bin/bash
CONF=/project/projectdirs/m779/yfeng1/imaginglss/dr2.conf.py

cp RANDOM.hdf5 QSO-RANDOM.hdf5

python ../scripts/imglss-query-tycho-veto.py QSO-RANDOM.hdf5 --conf $CONF
python ../scripts/imglss-query-tycho-veto.py QSO.hdf5 --conf $CONF
python ../scripts/imglss-query-completeness.py --conf $CONF \
        --use-tycho-veto=BOSS_DR9 QSO QSO.hdf5 QSO.hdf5
python ../scripts/imglss-query-completeness.py --conf $CONF \
        --use-tycho-veto=BOSS_DR9 QSO QSO.hdf5 QSO-RANDOM.hdf5
python '../scripts/imglss-export-text.py' --conf $CONF \
        --use-tycho-veto=BOSS_DR9 --bands r g z W1 W2 -- QSO.hdf5 QSO.txt
python '../scripts/imglss-export-text.py' --conf $CONF \
        --use-tycho-veto=BOSS_DR9 --bands r g z W1 W2 -- QSO-RANDOM.hdf5 QSO-RANDOM.txt

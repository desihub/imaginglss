#!/bin/bash
CONF=/project/projectdirs/m779/yfeng1/imaginglss/dr2.conf.py

cp RANDOM.hdf5 ELG-RANDOM.hdf5

python ../scripts/imglss-query-tycho-veto.py ELG-RANDOM.hdf5 --conf $CONF
python ../scripts/imglss-query-tycho-veto.py ELG.hdf5 --conf $CONF
python ../scripts/imglss-query-completeness.py --conf $CONF \
        --use-tycho-veto=BOSS_DR9 ELG ELG.hdf5 ELG.hdf5
python ../scripts/imglss-query-completeness.py --conf $CONF \
        --use-tycho-veto=BOSS_DR9 ELG ELG.hdf5 ELG-RANDOM.hdf5
python '../scripts/imglss-export-text.py' --conf $CONF \
        --use-tycho-veto=BOSS_DR9 --bands r g z W1 W2 -- ELG.hdf5 ELG.txt
python '../scripts/imglss-export-text.py' --conf $CONF \
        --use-tycho-veto=BOSS_DR9 --bands r g z W1 W2 -- ELG-RANDOM.hdf5 ELG-RANDOM.txt


python scripts/imglss-mpi-build-cache.py --conf=testdata/dr2-mini/dr2.conf.py || exit 1

python scripts/imglss-mpi-select-objects.py --conf=testdata/dr2-mini/dr2.conf.py ELG ELG.hdf5|| exit 1

python scripts/imglss-mpi-make-random.py --conf=testdata/dr2-mini/dr2.conf.py 10000 ELG-random.hdf5 || exit 1
python scripts/imglss-mpi-query-depth.py --conf=testdata/dr2-mini/dr2.conf.py ELG-random.hdf5 || exit 1

python scripts/imglss-query-tycho-veto.py --conf=testdata/dr2-mini/dr2.conf.py ELG.hdf5 || exit 1
python scripts/imglss-query-tycho-veto.py --conf=testdata/dr2-mini/dr2.conf.py ELG-random.hdf5 || exit 1

python scripts/imglss-query-completeness.py --conf=testdata/dr2-mini/dr2.conf.py --use-tycho-veto=BOSS_DR9 ELG ELG.hdf5 ELG-random.hdf5 || exit 1
python scripts/imglss-query-completeness.py --conf=testdata/dr2-mini/dr2.conf.py --use-tycho-veto=BOSS_DR9 ELG ELG.hdf5 ELG.hdf5 || exit 1

python scripts/imglss-export-text.py --conf=testdata/dr2-mini/dr2.conf.py --use-tycho-veto=BOSS_DR9 --bands g r z W1 -- ELG ELG.hdf5 ELG.txt
python scripts/imglss-export-text.py --conf=testdata/dr2-mini/dr2.conf.py --use-tycho-veto=BOSS_DR9 --bands g r z W1 -- ELG ELG-random.hdf5 ELG-random.txt
python scripts/imglss-naive-correlation.py --use-tycho-veto=BOSS_DR9 ELG.hdf5 ELG-random.hdf5 ELG-w.txt
python scripts/imglss-naive-crosscorrelation.py --use-tycho-veto=BOSS_DR9 ELG.hdf5 ELG-random.hdf5 ELG.hdf5 ELG-w.txt

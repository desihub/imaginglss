
python scripts/imglss-build-cache.py --conf=testdata/dr3-mini/dr3.conf.py || exit 1

python scripts/imglss-mpi-select-objects.py --conf=testdata/dr3-mini/dr3.conf.py LRG LRG.hdf5|| exit 1

python scripts/imglss-mpi-make-random.py --conf=testdata/dr3-mini/dr3.conf.py 10000 LRG-random.hdf5 || exit 1
python scripts/imglss-mpi-query-depth.py --conf=testdata/dr3-mini/dr3.conf.py LRG-random.hdf5 || exit 1

python scripts/imglss-query-tycho-veto.py --conf=testdata/dr3-mini/dr3.conf.py LRG.hdf5 || exit 1
python scripts/imglss-query-tycho-veto.py --conf=testdata/dr3-mini/dr3.conf.py LRG-random.hdf5 || exit 1

python scripts/imglss-query-completeness.py --conf=testdata/dr3-mini/dr3.conf.py --use-tycho-veto=BOSS_DR9 LRG LRG.hdf5 LRG-random.hdf5 || exit 1
python scripts/imglss-query-completeness.py --conf=testdata/dr3-mini/dr3.conf.py --use-tycho-veto=BOSS_DR9 LRG LRG.hdf5 LRG.hdf5 || exit 1

python scripts/imglss-export-text.py --conf=testdata/dr3-mini/dr3.conf.py --use-tycho-veto=BOSS_DR9 --bands g r z W1 -- LRG LRG.hdf5 LRG.txt
python scripts/imglss-export-text.py --conf=testdata/dr3-mini/dr3.conf.py --use-tycho-veto=BOSS_DR9 --bands g r z W1 -- LRG LRG-random.hdf5 LRG-random.txt
python scripts/imglss-naive-correlation.py --use-tycho-veto=BOSS_DR9 LRG.hdf5 LRG-random.hdf5 LRG-w.txt
python scripts/imglss-naive-crosscorrelation.py --use-tycho-veto=BOSS_DR9 LRG.hdf5 LRG-random.hdf5 LRG.hdf5 LRG-w.txt


python scripts/imglss-mpi-build-cache.py --conf=testdata/dr2-mini/dr2.conf.py || exit 1

python scripts/imglss-mpi-select-objects.py --conf=testdata/dr2-mini/dr2.conf.py QSO QSO.hdf5|| exit 1

python scripts/imglss-mpi-make-random.py --conf=testdata/dr2-mini/dr2.conf.py 10000 QSO-random.hdf5 || exit 1

python scripts/imglss-query-tycho-veto.py --conf=testdata/dr2-mini/dr2.conf.py QSO.hdf5 || exit 1
python scripts/imglss-query-tycho-veto.py --conf=testdata/dr2-mini/dr2.conf.py QSO-random.hdf5 || exit 1

python scripts/imglss-query-completeness.py --conf=testdata/dr2-mini/dr2.conf.py --use-tycho-veto=BOSS_DR9 QSO QSO.hdf5 QSO-random.hdf5 || exit 1
python scripts/imglss-query-completeness.py --conf=testdata/dr2-mini/dr2.conf.py --use-tycho-veto=BOSS_DR9 QSO QSO.hdf5 QSO.hdf5 || exit 1

python scripts/imglss-export-text.py --conf=testdata/dr2-mini/dr2.conf.py --use-tycho-veto=BOSS_DR9 --bands W1 -- QSO.hdf5 QSO.txt
python scripts/imglss-export-text.py --conf=testdata/dr2-mini/dr2.conf.py --use-tycho-veto=BOSS_DR9 --bands W1 -- QSO-random.hdf5 QSO-random.txt

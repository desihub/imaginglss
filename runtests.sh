
python scripts/imglss-mpi-build-cache.py --conf=testdata/dr2-mini/dr2.conf.py || exit 1
python scripts/imglss-mpi-select-objects.py --conf=testdata/dr2-mini/dr2.conf.py QSO QSO.txt || exit 1
head QSO.txt.NOISES
head QSO.txt.FLUXES
head QSO.txt.CONFIDENCE

python scripts/imglss-mpi-select-objects.py --conf=testdata/dr2-mini/dr2.conf.py QSO QSO.hdf5 || exit 1
python scripts/imglss-mpi-make-random.py --conf=testdata/dr2-mini/dr2.conf.py 1000 QSO-random.txt || exit 1
head QSO-random.txt.NOISES
python scripts/imglss-mpi-make-random.py --conf=testdata/dr2-mini/dr2.conf.py 1000 QSO-random.hdf5 || exit 1

python scripts/imglss-query-completeness.py --conf=testdata/dr2-mini/dr2.conf.py QSO QSO-random.txt QSO.txt QSO-random.txt || exit 1
head QSO-random.txt.FC
python scripts/imglss-query-tycho-veto.py --conf=testdata/dr2-mini/dr2.conf.py QSO-random.txt QSO-random.txt|| exit 1
head QSO-random.txt.TYCHOVETO

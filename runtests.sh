
python scripts/build_cache.py --conf=testdata/dr2-mini/dr2.conf.py || exit 1
python scripts/select_objs.py --conf=testdata/dr2-mini/dr2.conf.py --with-tycho=DECAM_QSO QSO QSO.txt || exit 1
head QSO.txt
python scripts/select_objs.py --conf=testdata/dr2-mini/dr2.conf.py --with-tycho=DECAM_QSO QSO QSO.hdf5 || exit 1
python scripts/make_random.py --conf=testdata/dr2-mini/dr2.conf.py --with-tycho=DECAM_QSO 1000 QSO QSO-random.txt || exit 1
head QSO-random.txt
python scripts/make_random.py --conf=testdata/dr2-mini/dr2.conf.py --with-tycho=DECAM_QSO 1000 QSO QSO-random.hdf5 || exit 1

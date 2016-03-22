
python scripts/build_cache.py --conf=testdata/dr2-mini/dr2.conf.py || exit 1
python scripts/select_objs.py --conf=testdata/dr2-mini/dr2.conf.py --with-tycho=DECAM_QSO QSO QSO.txt || exit 1
head QSO.txt.NOISES
head QSO.txt.FLUXES
python scripts/select_objs.py --conf=testdata/dr2-mini/dr2.conf.py --with-tycho=DECAM_QSO QSO QSO.hdf5 || exit 1
python scripts/make_random.py --conf=testdata/dr2-mini/dr2.conf.py --with-tycho=DECAM_QSO 1000 QSO-random.txt || exit 1
head QSO-random.txt.NOISES
python scripts/make_random.py --conf=testdata/dr2-mini/dr2.conf.py --with-tycho=DECAM_QSO 1000 QSO-random.hdf5 || exit 1

python scripts/query_completeness.py --conf=testdata/dr2-mini/dr2.conf.py QSO QSO-random.txt QSO.txt QSO-random.txt || ext 1
python scripts/query_veto.py --conf=testdata/dr2-mini/dr2.conf.py QSO-random.txt QSO-random.txt|| ext 1

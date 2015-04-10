# bash

if ! python -c 'import numpydoc'; then easy_install --user numpydoc; fi
if ! python -c 'import sphinx'; then easy_install --user sphinx; fi

sphinx-apidoc -M -e -f -o . ../model ../model/utils/sharedmem

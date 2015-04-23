# bash

if ! python -c 'import numpydoc'; then easy_install --user numpydoc; fi
if ! python -c 'import sphinx'; then easy_install --user sphinx; fi

sphinx-apidoc -H "API Reference: model" -M -e -f -o . ../model ../model/utils/sharedmem
sphinx-apidoc -H "API Reference: analysis" -M -e -f -o . ../analysis

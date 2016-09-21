from setuptools import setup, find_packages
from glob import glob
setup(name="imaginglss", version="0.1.0rc0",
      author="Yu Feng, Martin White, Ellie Kitanidis",
      maintainer="Yu Feng",
      maintainer_email="rainwoodman@gmail.com",
      description="Generating large scale structure catalogues from DECALS data",
      url="http://github.com/desihub/imaginglss",
      zip_safe=False,
      package_dir = {'imaginglss': 'imaginglss'},
      include_package_data=False,
      packages = find_packages(),
      scripts = [
            'scripts/imglss-build-cache.py',
            'scripts/imglss-export-text.py',
            'scripts/imglss-mpi-make-random.py',
            'scripts/imglss-mpi-query-depth.py',
            'scripts/imglss-mpi-select-objects.py',
            'scripts/imglss-query-completeness.py',
            'scripts/imglss-query-tycho-veto.py',
            'scripts/imglss-naive-correlation.py',
            'scripts/imglss-naive-crosscorrelation.py',
      ],
      install_requires=['numpy'],
      requires=['bigfile'],
)


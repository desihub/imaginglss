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
      install_requires=['numpy']
)


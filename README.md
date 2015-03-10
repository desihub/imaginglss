# mydesi

This is a data model of the DECAM imaging survey data 
and potentially also for other desi targeting imaging survey, 
for Large scale structure analysis.

The modelled data is a product from the Tractor, the source identification
software used in by DESI.

We have support on the Early Data Release (EDR) data, and DR1 support is on the way.

See 
https://desi.lbl.gov/trac/wiki/DecamLegacy/EDRfiles

and 
https://desi.lbl.gov/trac/wiki/DecamLegacy/DR1

To use this imaing data model on edison:

1. set environments
  
  ```
    export DECALS_IMAGING=/global/project/projectdirs/cosmo/work/decam/release/edr/
    export DECALS_CACHE=$GSCRATCH/desicache
  ```
2. install astropy or fitsio (prefered)
  ```
    easy_install --user astropy
  ```

  Alternatively, 
  ```
    easy_install --user fitsio
  ```

3. try if usecase.py runs.
  ```
    python usecase.py
  ```

4. start hacking with
  ```
    from model.datarelease import DataRelease
    dr = DataRelease()
  ```

Due to internal caching of the catalogue, the running time is faster after first time.

Documentation of the modelling can be seen in model/README.md

We are working on a more comprehensive documentation generated from
source code doc strings.

Examples:
  ```
  edison_example.sh
  mpirandom.py
  usecase.py
  select_elgs.py
  footprint_mangle.py
  build_sfdmap.py
  ```


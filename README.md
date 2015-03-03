# mydesi

To use the data model on edison:

1. set environments:


    export DECALS_IMAGING=/global/project/projectdirs/cosmo/work/decam/release/edr/
    export DECALS_CACHE=$HOME/desicache

2. somehow install astropy.

3. try if usecase.py runs.
    

Describe the derived products that worth saving.

* depths in r, u, z bands:

   reading off of depth-\*-[ruz].fits.gz for object
   in tractor-edr. 

   It takes a few minitues for the 3 million objects in tractor-edr.fits

   In the future, the catalogue comes in tractor directory
   as one file per brick. We need a smarter way of handling 
   this.

* bare bone catalogue files.
  
  RA DEC and r, u, z nanomags.

* sky brightness:

   Where shall this come from?
   is this also in bricks?

* extinction:
   
   extinction in r, u, z (and maybe more) as a function of RA DEC;
   is this also in bricks?
   is it on Nersc

* TACO star catalog
   is it on nersc?

* Details on CFHTLS can be found at:
http://www.cfht.hawaii.edu/Science/CFHTLS/

there was also a series of nice papers by the "CARS" collaboration (The CFHTLS-Archive-Research Survey) which I remember being quite readable.  They also made a lot of processed data available.  The following looks useful:

http://adsabs.harvard.edu/abs/2009A%26A...493.1197E

and might give us some ideas for figures/tests/cross-checks we can do.

* DESI software on Edison: /project/projectdirs/desi/software/edison/
  dust directory has dust

* depth as function of RA dEC in mollleview??

* compare SDSS flux with DECAM flux.
  fit a line. 

"""
Python code to read and write "filehander files", which are
simply directories with files in them for each "keyword".


"""
from __future__ import print_function

__author__ = "Martin White and Yu Feng"
__version__ = "1.0"
__email__  = "mwhite@berkeley.edu yfeng1@berkeley.edu"


import numpy as N
import glob
import re
import os
import os.path

def read_file(fname,keys=None):
    """
    Reads a specificed list of "keys" in a "file" of name fname.

    Parameters
    ----------
    fname : string
        location to look for the data columns
    keys : list or None
        list of the keys; all if keys==None.  

    Returns
    -------
    data : dict
        a dictionary of NumPy arrays.

    Notes
    -----
    Does minimal checking, assuming you know what you're doing.
    The on-disk representation is always in '<'. We swap native endianness, on the fly if it
    differ from the native endianness.

    """
    # Get a list of all of the "fields" in the "file".
    flist = glob.glob(fname+"/*[fi][48]")
    # and start filling in my dictionary.
    ret = {}

    if keys is None:
        keys = [parse_filename(os.path.basename(fn))[0] for fn in flist]

    for fn in flist:
        key, dtype = parse_filename(os.path.basename(fn))
        # and add it to the dictionary
        if key not in keys: continue
        order = dtype.str[0]
        with open(fn, "r") as ff:
            d = N.fromfile(ff, dtype=dtype)
        if order != '<':
            d.byteswap(True)
        ret[key] = d

    # Now check we got everything.
    for key in keys:
        if key in ret.keys(): continue
        raise RuntimeError("Unable to find "+key+" in "+fname)

    return(ret)
    #

def format_filename(key, data):
    """ Generate a filename from base name for data 

        Parameters
        ----------
        data : array_like
            Must be 1d or 2d array. For 2d array, the shape
            of the last dimension is used as the size of the vector.

        key: string
            key
        
        Returns
        -------
        filename : string
            key.vectorsize.datatype
            
    """
    assert len(data.shape) <= 2
    if len(data.shape) == 2:
       vectorsize = data.shape[-1] 
    elif len(data.shape) == 1:
       vectorsize = 0

    suffix={}
    suffix['int32']='i4'
    suffix['int64']='i8'
    suffix['float32']='f4'
    suffix['float64']='f8'

    return "%s.%d.%s" % (key,
        vectorsize,
        suffix[data.dtype.name])

def parse_filename(filename):
    """ 
        parse a single file name

        Parameters
        ----------
        filename : string
            filename
        
        Returns
        -------
        key: string
            key of the 
        dtype: dtype
            dtype

    """
    mm = re.search(r"(\w*)\.([0-9]+)\.([fi][48])", filename)

    if mm==None:
        raise RuntimeError("Unable to parse file name `%s`." % filename)
    else:
        key = mm.group(1)
        vectorsize = int(mm.group(2))
        # force the endianness
        objt= N.dtype(mm.group(3))

    if vectorsize > 0:
        objt = N.dtype((objt, vectorsize))
        
    return key, objt

def write_file(fname, data, mode='w'):
    """
    Writes the dictionary, data, which is meant to contain only
    NumPy arrays, to a "file" of name fname. The file is always
    written with '<' endian.

    Parameters
    ----------
    fname : string
        location where the files will be written
    data : dict alike
        key and array pairs. Each item must be a numerical numpy array
    mode : string 
        'w' for writing or 'a' for appending.
    
    Notes
    -----
    Does minimal checking, assuming you know what you're doing.


    """
    # Put in a little object type converter.
    if not os.path.exists(fname):
        os.mkdir(fname)
    for key in data.keys():
        filename = os.path.join(
                fname, format_filename(key, data[key]))
        d = data[key]
        order = d.dtype.str[0]
        if order != '<':
            d.byteswap(True)
        with open(filename, mode) as ff:
            d.tofile(ff)
        if order != '<':
            d.byteswap(True)
    #

def test():
    a = {
        'a': N.arange(10),
        'b': N.arange(10, dtype='>f4'),
    }
    write_file('filehandler-test', a, mode='w')
    write_file('filehandler-test', a, mode='a')
    data = read_file('filehandler-test')
    assert len(data['a']) == 20
    assert len(data['b']) == 20
    print(data)

if __name__ == '__main__':
    test()

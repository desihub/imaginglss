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
class MissingColumn(IOError):
    pass
class BadFileName(IOError):
    pass

def read(fname, keys=None, offset=0, count=None):
    """
    Reads a specificed list of "keys" in a "file" of name fname.

    Parameters
    ----------
    fname : string
        location to look for the data columns
    keys : list or None
        list of the keys; all if keys==None.  
    offset : int
        offset to start reading, in unit of items
    count  : int
        total number to read. None for read to the end of the file

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
    # and start filling in my dictionary.
    ret = {}

    columns = list(fname)
    if keys is None:
        keys = columns.keys()

    for key in keys:
        dtype = columns[key]
        # and add it to the dictionary
        order = dtype.base.str[0]
        fn = format_filename(key, dtype)
        with open(os.path.join(fname, fn), "r") as ff:
            ff.seek(dtype.itemsize * offset, 0)
            if count is None:
                count = -1
            d = N.fromfile(ff, dtype=dtype, count=count)
        if order != '<':
            d.byteswap(True)
        ret[key] = d

    # Now check we got everything.
    for key in keys:
        if key in ret.keys(): continue
        raise MissingColumn("Unable to find "+key+" in "+fname)

    return(ret)
    #

def list(fname):
    # Get a list of all of the "fields" in the "file".
    flist = glob.glob(fname+"/*.*.*")
    # and start filling in my dictionary.
    ret = {}

    for fn in flist:
        key, dtype = parse_filename(os.path.basename(fn))
        ret[key] = dtype
    return ret

def size(fname, key):
    """ Number of items in this column """

    columns = list(fname)
    if key not in columns:
        raise MissingColumn("Unable to find "+key+" in "+fname)

    dtype = columns[key]
    size = os.path.getsize(os.path.join(fname, format_filename(key, dtype)))
    assert size % dtype.itemsize == 0
    return size // dtype.itemsize

def build_dtype(data):
    if len(data.shape) >= 2:
        return N.dtype((data.dtype, data.shape[1:]))
    elif len(data.shape) == 1:
        return data.dtype
    else:
        raise TypeError("data is not at least 1d")

def format_filename(key, data_or_dtype):
    """ Generate a filename from base name for data 

        Parameters
        ----------
        data_or_dtype : array_like or dtype
            Must be 1d or 2d array. For 2d array, the shape
            of the last dimension is used as the size of the vector.
            Or a dtype with shape

        key: string
            key
        
        Returns
        -------
        filename : string
            key.vectorsize.datatype
            
    """
    if not isinstance(data_or_dtype, N.dtype):
        dtype = build_dtype(data_or_dtype)
    else:
        dtype = data_or_dtype

    if len(dtype.shape) > 0:
        vectorsize = 'x'.join(['%d' % x for x in dtype.shape[:]])
    else:
        vectorsize = '0'

    dtype = dtype.base

    suffix = dtype.str[1:]

    return "%s.%s.%s" % (key,
        vectorsize,
        suffix)

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
    mm = re.search(r"(\w*)\.([0-9x]+)\.([a-zA-Z][0-9]*)", filename)

    if mm==None:
        raise BadFileName("Unable to parse file name `%s`." % filename)
    else:
        key = mm.group(1)
        vectorsize = [int(x) for x in mm.group(2).split('x')]
        # force the endianness
        objt= N.dtype(mm.group(3))

    if not (len(vectorsize) == 1 and vectorsize[0] == 0):
        objt = N.dtype((objt, tuple(vectorsize)))
    return key, objt

def write(fname, data, mode='w', offset=None):
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
        'w' for writing or 'a' for appending, or 'r+' for writing at
        offset.
    offset : int or None
        if the mode is 'r+',
        offset in a file for the write operation, in units of items.
        None to append to the end.
 
    Notes
    -----
    Does minimal checking, assuming you know what you're doing.


    """
    if offset is not None:
        if mode != 'r+':
            raise ValueError("To use an offset, mode must be 'r+'")

    # Put in a little object type converter.
    if not os.path.exists(fname):
        os.makedirs(fname)

    if isinstance(data, N.ndarray):
        keys = data.dtype.names
    else:
        keys = data.keys()
        
    for key in keys:
        filename = os.path.join(
                fname, format_filename(key, data[key]))
        d = data[key]
        order = d.dtype.base.str[0]
        if order != '<':
            d.byteswap(True)
        dtype = build_dtype(d)
        with open(filename, mode) as ff:
            if offset is not None:
                ff.seek(offset * dtype.itemsize, 0)
            d.tofile(ff)

        if order != '<':
            d.byteswap(True)
    #

def test():
    a = {
        'a': N.arange(10),
        'b': N.arange(10, dtype='>f4'),
        'c': N.arange(10, dtype='>i4').reshape(2, 5),
    }
    write('filehandler-test', a, mode='w')
    write('filehandler-test', a, mode='a')
    write('filehandler-test', a, mode='r+', offset=20)
    data = read('filehandler-test')
    assert len(data['a']) == 30
    assert len(data['b']) == 30
    assert len(data['c']) == 22
    print(data)

if __name__ == '__main__':
    test()

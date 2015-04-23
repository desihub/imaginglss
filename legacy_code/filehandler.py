#
# Python code to read and write "filehander files", which are
# simply directories with files in them for each "keyword".
#
from __future__ import print_function

__author__ = "Martin White"
__version__ = "1.0"
__email__  = "mwhite@berkeley.edu"


import numpy as N
import glob
import re
import os


def read_file(fname,keys=None):
    """
    read_file(fname,keys=None):
    Reads a specificed list of "keys" in a "file" of name fname,
    or all of the keys if keys==None.  Returns a dictionary of
    NumPy arrays.
    Does minimal checking, assuming you know what you're doing.
    """
    # Get a list of all of the "fields" in the "file".
    flist = glob.glob(fname+"/*[fi][48]")
    # and start filling in my dictionary.
    ret = {}
    for fn in flist:
        # Get the field name and type.
        mm = re.search(fname+"/"+r"(\w*)\.([fi][48])",fn)
        if mm==None:
            raise RuntimeError,"Unable to parse file "+fn
        else:
            key = mm.group(1)
            objt= mm.group(2)
        objt = '<' + objt
        # and add it to the dictionary
        if keys==None:	# Need to do it this way since can't iterate None.
            ff  = open(fn,"r")
            nobj= N.fromfile(ff,count=1,dtype='<i8')
            ret[key]=N.fromfile(ff,count=nobj,dtype=objt)
            ff.close()
        elif key in keys:
            ff  = open(fn,"r")
            nobj= N.fromfile(ff,count=1,dtype='<i8')
            ret[key]=N.fromfile(ff,count=nobj,dtype=objt)
            ff.close()
    # Now check we got everything.
    if keys!=None:
        for key in keys:
            if key not in ret.keys():
                raise RuntimeError,"Unable to find "+key+" in "+fname
    return(ret)
    #



def write_file(fname,data):
    """
    write_file(fname,data):
    Writes the dictionary, data, which is meant to contain only
    NumPy arrays, to a "file" of name fname.
    Does minimal checking, assuming you know what you're doing.
    """
    # Put in a little object type converter.
    suffix={}
    suffix['int32']='i4'
    suffix['int64']='i8'
    suffix['float32']='f4'
    suffix['float64']='f8'
    if not os.path.exists(fname):
        os.mkdir(fname)
    for key in data.keys():
        dt  = suffix[data[key].dtype.name]
        ff  = open(fname+"/"+key+"."+dt,"w")
        nobj= N.array( data[key].size, dtype='<i8' )
        nobj.tofile(ff)
        data[key].astype('<'+dt).tofile(ff)
        ff.close()
    #

import numpy
import os
class writer(object):
    def __init__(self, filename):
        if filename.endswith('.txt'):
            self.write_output = write_text
        elif filename.endswith('.hdf5'):
            import h5py
            self.write_output = write_hdf5
        else :
            raise ValueError('Output must be either txt or hdf5')
        self.filename = filename

    def write(self, data, metadict):
        self.write_output(self.filename, data, metadict)

    def __str__(self):
        return os.path.abspath(self.filename)

def write_text(output, CANDIDATES, metadict):
    names = CANDIDATES.dtype.names
    def format_name(name, dtype):
        subdtype = dtype[name]
        if subdtype.shape is not None and len(subdtype.shape):
            return '%s[%s]' % (name, subdtype.shape)
        else:
            return name
    def format_var(var, fmt, subdtype):
        if subdtype.shape is not None and len(subdtype.shape):
            return ' '.join([fmt % var[i] for i in range(subdtype.shape[0])])
        else:
            return fmt % var

    with open(output, "w") as ff:
        ff.write('\n'.join(['# %s = %s' % (str(key), str(metadict[key])) for key in sorted(metadict.keys())]))
        ff.write("\n# %s\n" % ' '.join([format_name(name, CANDIDATES.dtype) for name in names]))
        for row in CANDIDATES:
            for name in names:
                ff.write(format_var(row[name], '%15.10f ', CANDIDATES.dtype[name]))
            ff.write('\n')

def write_hdf5(output, CANDIDATES, metadict):
    import h5py
    with h5py.File(output, 'w') as ff:
        ds = ff.create_dataset('CANDIDATES', data=CANDIDATES)
        for key, v in metadict.items():
            if isinstance(v, object):
                ds.attrs[key] = str(v)
            else:
                ds.attrs[key] = v


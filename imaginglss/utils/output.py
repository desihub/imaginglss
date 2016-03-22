import numpy
import os
class writer(object):
    def __init__(self, filename):
        if filename.endswith('.txt'):
            self.write_output = write_text
            self.read_output = read_text
        elif filename.endswith('.hdf5'):
            import h5py
            self.write_output = write_hdf5
            self.read_output = read_hdf5
        else :
            raise ValueError('Output must be either txt or hdf5')
        self.filename = filename

    def write(self, data, metadict, extension):
        self.write_output(self.filename, data, metadict, extension)

    def read(self, extension):
        return self.read_output(self.filename, extension)

    def __str__(self):
        return os.path.abspath(self.filename)

def read_text(output, extension):
    dtypes = {
    'FLUXES' : numpy.dtype([
        ('RA', 'f8'),
        ('DEC', 'f8'),
        ('PHOTO_Z', 'f8'),
        ('DECAM_INTRINSIC_FLUX', ('f8', 6))]),

    'NOISES' : numpy.dtype([
        ('RA', 'f8'),
        ('DEC', 'f8'),
        ('DECAM_INTRINSIC_NOISE_LEVEL', ('f8', 6)),
        ])
    }
    dtype = dtypes[extension]

    data = numpy.loadtxt(output + '.' + extension, dtype='f8')
    data = data.view(dtype=dtype).reshape(-1)
    return data

def write_text(output, CANDIDATES, metadict, extension):
    names = CANDIDATES.dtype.names
    fmts = {'f': '%15.10f ', 'd': '%15.10f ', 'u': '0X%013X ', 'i' : '%15d ', 'b' : '%01d '}

    def format_name(name, dtype):
        subdtype = dtype[name]
        if subdtype.shape is not None and len(subdtype.shape):
            return '%s[%s]' % (name, subdtype.shape)
        else:
            return name
    def format_var(var, subdtype):
        if subdtype.shape is not None and len(subdtype.shape):
            fmt = fmts[subdtype.base.kind]
            return ' '.join([fmt % var[i] for i in range(subdtype.shape[0])])
        else:
            fmt = fmts[subdtype.kind]
            return fmt % var

    with open(output + '.' + extension, "w") as ff:
        ff.write('\n'.join(['# %s = %s' % (str(key), str(metadict[key])) for key in sorted(metadict.keys())]))
        ff.write("\n# %s\n" % ' '.join([format_name(name, CANDIDATES.dtype) for name in names]))
        for row in CANDIDATES:
            for name in names:
                char = CANDIDATES.dtype[name].char
                ff.write(format_var(row[name], CANDIDATES.dtype[name]))
            ff.write('\n')

def read_hdf5(output, extension):
    import h5py
    with h5py.File(output) as ff:
        return ff[extension][:]

def write_hdf5(output, CANDIDATES, metadict, extension):
    import h5py
    with h5py.File(output) as ff:
        if extension in ff:
            del ff[extension]
        ds = ff.create_dataset(extension, data=CANDIDATES)
        for key, v in metadict.items():
            if isinstance(v, object):
                ds.attrs[key] = str(v)
            else:
                ds.attrs[key] = v


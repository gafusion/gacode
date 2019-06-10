from gacode import expro
import numpy
import sys
if sys.version_info < (3, 0):
    def b2s(obj):
        return obj
else:
    def b2s(obj):
        if isinstance(obj, bytes):
            return obj.decode("utf-8")
        import numpy
        if isinstance(obj, numpy.ndarray) and obj.dtype.name.startswith('bytes'):
            return numpy.reshape(numpy.array(list(map(b2s, obj.flat))), obj.shape)
        else:
            return obj

class ExproData(dict):

    def __init__(self, filename, input_profiles_compatibility_mode=True):
        expro.expro_read(filename)
        for item in dir(expro):
            if item.startswith('expro_') and item.endswith('_str'):
                try:
                    self[item[len('expro_'):-len('str_')]] = getattr(expro, item[:-len('str_')])
                except Exception as _excp:
                    print('WARNING: issue reading %s: %s' % (item, repr(_excp)))

        if input_profiles_compatibility_mode:

            self['Te'] = self['te']
            del self['te']
            self['Ti'] = self['ti']
            del self['ti']
            for quantity in ['ni', 'Ti', 'vpol', 'vtor']:
                if quantity in self:
                    for k in range(self['n_ion']):
                        self[quantity + '_%d' % (k + 1)] = self[quantity][k]
                    del self[quantity]
            self['IONS'] = []
            for k in range(self['n_ion']):
                self['IONS'].append([b2s(self['name'][k]).strip(), self['z'][k], self['mass'][k], b2s(self['type'][k]).strip()])
            del self['z']
            del self['mass']
            del self['ze']
            del self['masse']
            print(self['IONS'])

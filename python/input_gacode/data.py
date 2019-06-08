from gacode import expro
import numpy

class ExproData(dict):

    def __init__(self, filename, input_profiles_compatibility_mode=True):
        expro.expro_read(filename)
        for item in dir(expro):
            if item.startswith('expro_') and item.endswith('_str'):
                try:
                    self[item[len('expro_'):-len('str_')]] = getattr(expro, item[:-len('str_')])
                except Exception as _excp:
                    print('WARNING: issue reading %s: %s' % (item, repr(_excp)))
        for item in ['name', 'type']:
            self[item] = map(lambda x: ''.join(x).strip().replace('[', '').replace(']', ''),
                             numpy.split(numpy.array(list(self[item].tolist())), self['n_ion']))

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
                self['IONS'].append([self['name'][k], self['z'][k], self['mass'][k], self['type'][k]])
            del self['z']
            del self['mass']
            del self['ze']
            del self['masse']
            print(self['IONS'])

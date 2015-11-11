#!/opt/local/bin/python

with open('EXPRO_interface.f90','r') as f:
    fl = f.readlines()
var_names = []
tags = {}
dims = {}
for l in fl:
    ll = l.lower().strip().split('!')[0]
    if not ll:
        continue
    if '::' in ll:
        dim = 0
        if 'dimension' in ll:
            dim = ll.split('::')[0].count(':')
        for var0 in ll.split('::')[1].split(','):
            var = var0.strip()
            if '=' in var:
                var_name,val = map(str.strip,var.split('='))
            else:
                var_name = var.strip()
                val = 'not set'
            if var_name.endswith('_tag') and var_name[:-4] in var_names:
                try:
                    tags[var_name[:-4]] = eval(val)
                except:
                    tags[var_name[:-4]] = val
            else:
                var_names.append(var_name)
                dims[var_name] = dim
quants_per_ion = ['ti','vpol','ni','vtor']
print var_names
print tags
units = {
            'ti': 'keV',
            'ni': '10^19/m^3',
            'vpol': 'm/s',
            'vtor': 'm/s'
        }
for k,v in tags.items():
    tag,unit = v.split('(',2)[0:2]
    unit = ')'.join(unit.split(')')[:-1])
    if k[len('expro_'):]!=tag:
        print k,tag
        continue
    if tag in units:
        continue
    units[tag] = unit 
    if '(/' in v and '(/m' not in v:
        print k
print units
print dims
with open('expro_vars.txt','w') as f:
    for k in var_names:
        f.write('----\n')
        f.write('Variable: %s\n'%k.replace('expro_','EXPRO_'))
        f.write('Description: \n')
        f.write('Units: %s\n'%units.get(k.replace('expro_',''),'?').replace('mw','MW').replace('kev','keV'))
        f.write('Dimensions: %d\n'%dims[k])
        if dims[k]>0:
            f.write('Dimension 1: rho\n')
            if dims[k]>1:
                f.write('Dimension 2: ions\n')
                if dims[k]>2:
                    for d in range(3,dims[k]+1):
                        f.write('Dimension %d:?\n'%d)

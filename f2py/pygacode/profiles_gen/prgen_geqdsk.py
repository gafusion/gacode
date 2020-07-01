from scipy import interpolate, integrate
import numpy as np
import fortranformat
import os

def splitter(inv, step=16):
    value = []
    for k in range(int(len(inv) / step)):
        value.append(inv[step * k:step * (k + 1)])
    return value

def merge(inv):
    return ''.join(inv)

def prgen_geqdsk(filename):

    slf = {}
    with open(filename, 'r') as f:
        EQDSK = f.read().splitlines()

    slf['CASE'] = np.array(splitter(EQDSK[0][0:48],8))

    try:
        tmp = list([_f for _f in EQDSK[0][48:].split(' ') if _f])
        [IDUM, slf['NW'],slf['NH']] = list(map(int,tmp[:3]))
    except ValueError: # Can happen if no space between numbers, such as 10231023
        IDUM = int(EQDSK[0][48:52])
        slf['NW'] = int(EQDSK[0][52:56])
        slf['NH'] = int(EQDSK[0][56:60])
        tmp = []
        printd('IDUM,NW,NH',IDUM,slf['NW'],slf['NH'],topic='OMFITgeqdsk.load')

    if len(tmp) > 3:
        slf['EXTRA_HEADER'] = EQDSK[0][49+len(re.findall('%d +%d +%d ' % (IDUM,slf['NW'],slf['NH']),EQDSK[0][49:])[0])+2:]

    print('INFO: (prgen_geqdsk) Imported {:d}x{:d} gEQDSK'.format(slf['NW'],slf['NH']))

    offset = 1

    # now, the next 20 numbers (5 per row)
    [slf['RDIM'],   slf['ZDIM'],  slf['RCENTR'], slf['RLEFT'], slf['ZMID'], 
     slf['RMAXIS'], slf['ZMAXIS'],slf['SIMAG'],  slf['SIBRY'], slf['BCENTR'],
     slf['CURRENT'],slf['SIMAG'], XDUM,          slf['RMAXIS'],XDUM,
     slf['ZMAXIS'], XDUM,         slf['SIBRY'],  XDUM,         XDUM] \
     = list(map(eval,splitter(merge(EQDSK[offset:offset+4]))))

    offset = offset + 4

    # now I have to read NW elements
    nlNW = int(np.ceil(slf['NW']/5.0))

    slf['FPOL'] = np.array(list(map(float,splitter(merge(EQDSK[offset:offset+nlNW])))))
    offset = offset + nlNW
    slf['PRES'] = np.array(list(map(float,splitter(merge(EQDSK[offset:offset+nlNW])))))
    offset = offset + nlNW
    slf['FFPRIM'] = np.array(list(map(float,splitter(merge(EQDSK[offset:offset+nlNW])))))
    offset = offset + nlNW
    slf['PPRIME'] = np.array(list(map(float,splitter(merge(EQDSK[offset:offset+nlNW])))))
    offset = offset + nlNW

    try:
        # official gEQDSK file format saves PSIRZ as a single flat array of size rowsXcols
        nlNWNH = int(np.ceil(slf['NW']*slf['NH'] / 5.))
        slf['PSIRZ'] = np.reshape(np.fromiter(splitter(''.join(EQDSK[offset:offset+nlNWNH])),dtype=np.float),(slf['NH'],slf['NW']))
        offset = offset+nlNWNH
    except ValueError:
        # sometimes gEQDSK files save row by row of the PSIRZ grid (eg. FIESTA code)
        nlNWNH = slf['NH']*nlNW
        slf['PSIRZ'] = np.reshape(np.fromiter(splitter(''.join(EQDSK[offset:offset+nlNWNH])),dtype=np.float),(slf['NH'],slf['NW']))
        offset = offset+nlNWNH

    slf['QPSI'] = np.array(list(map(float,splitter(merge(EQDSK[offset:offset+nlNW])))))
    offset = offset+nlNW

    # now vacuum vessel and limiters
    [slf['NBBBS'],slf['LIMITR']] = list(map(int,[_f for _f in EQDSK[offset:offset + 1][0].split(' ') if _f]))
    offset = offset+1

    nlNBBBS = int(np.ceil(slf['NBBBS'] * 2 / 5.))
    slf['RBBBS'] = np.array(list(map(float,splitter(merge(EQDSK[offset:offset+nlNBBBS]))))[0::2])
    slf['ZBBBS'] = np.array(list(map(float,splitter(merge(EQDSK[offset:offset+nlNBBBS]))))[1::2])
    offset = offset+max(nlNBBBS,1)

    try:
        # this try/except is to handle some gEQDSK files written by older versions of ONETWO
        nlLIMITR = int(np.ceil(slf['LIMITR'] * 2 / 5.))
        slf['RLIM'] = np.array(list(map(float, splitter(merge(EQDSK[offset:offset + nlLIMITR]))))[0::2])
        slf['ZLIM'] = np.array(list(map(float, splitter(merge(EQDSK[offset:offset + nlLIMITR]))))[1::2])
        offset = offset + nlLIMITR
    except ValueError:
       print('ERROR: (prgen_geqdsk) EFIT Limiter error')
       
    try:
        [slf['KVTOR'], slf['RVTOR'], slf['NMASS']] = list(map(float,[_f for _f in EQDSK[offset:offset + 1][0].split(' ') if _f]))
        offset = offset + 1

        if slf['KVTOR'] > 0:
            slf['PRESSW'] = np.array(list(map(float, splitter(merge(EQDSK[offset:offset + nlNW])))))
            offset = offset + nlNW
            slf['PWPRIM'] = np.array(list(map(float, splitter(merge(EQDSK[offset:offset + nlNW])))))
            offset = offset + nlNW

        if slf['NMASS'] > 0:
            slf['DMION'] = np.array(list(map(float, splitter(merge(EQDSK[offset:offset + nlNW])))))
            offset = offset + nlNW

        slf['RHOVN'] = np.array(list(map(float, splitter(merge(EQDSK[offset:offset + nlNW])))))
        offset = offset + nlNW
    except Exception:
        pass

    # fix some gEQDSK files that do not fill PRES info (eg. EAST)
    if not np.sum(slf['PRES']):
        pres = integrate.cumtrapz(slf['PPRIME'], np.linspace(slf['SIMAG'], slf['SIBRY'], len(slf['PPRIME'])), initial=0)
        slf['PRES'] = pres - pres[-1]

    slf['PSI'] = np.linspace(slf['SIMAG'],slf['SIBRY'],len(slf['PRES']))
    
    return slf


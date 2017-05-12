#!/usr/bin/env python

def ddb_float(float_in):
    '''
    Convert float to Decimal compatible with DynamoDB format

    :param float_in: input float

    :return: float in Decimal format
    '''
    if not isinstance(float_in,basestring):
        float_in=repr(float(float_in))
    from boto3.dynamodb.types import DYNAMODB_CONTEXT
    return DYNAMODB_CONTEXT.create_decimal(float_in)

def _data_2_message(payload):
    import numpy
    import re
    def formatter(data):
        data=('%g'%data)
        if '.' in data:
            data=data.lstrip('0')
        data=re.sub('e\-0+','e-',data)
        data=re.sub('e\+0+','e+',data)
        return data
    def compress(tmps):
        return tmps
        tmpsc=[]
        c=-1
        kold=tmps[0]
        for k in tmps+[tmps[-1]+' ']:
            if kold==k:
                c+=1
            else:
                if c==0:
                    tmpsc.append(kold)
                else:
                    tmpsc.append(str(c+1)+'*'+kold)
                kold=k
                c=0
        return tmpsc

    message=[]
    for what in payload.keys():
        data=payload[what]
        if isinstance(data,bool):
            tp='b'
            data=str(int(data))
        elif isinstance(data,(list,tuple,numpy.ndarray)):
            tp='a'
            data=re.sub(' ','','['+','.join(compress(map(formatter,numpy.atleast_1d(data).flatten().tolist())))+']' )
        elif numpy.array(data).dtype.kind=='i':
            tp='i'
            data=str(data)
        elif numpy.array(data).dtype.kind=='f':
            tp='f'
            data=formatter(data)
        elif isinstance(data,basestring):
            tp='s'
            data=repr(data.strip())
        elif data is None:
            tp='s'
            data=''
        else:
            raise(Exception('%s objects of type %s are not supported'%(what,type(data))))
        message.append(tp+'@'+what+'='+data)

    return '|'.join(message)

def harvest_send(payload, table='test_harvest', host=None, port=None, verbose=None, tag=None, protocol=None):
    '''
    Function to send data to the harvesting server

    :param payload: dictionary

    :param table: table where to put the data

    :param host: harvesting server address
    If None take value from `HARVEST_HOST` environemental variable, or use default `gadb-harvest.ddns.net` if not set.

    :param port: port the harvesting server is listening on.
    If None take value from `HARVEST_PORT` environemental variable, or use default `0` if not set.

    :param verbose: print harvest message to screen
    If None take value from `HARVEST_VERBOSE` environemental variable, or use default `False` if not set.

    :param tag: tag entry
    If None take value from `HARVEST_TAG` environemental variable, or use default `Null` if not set.

    :param protocol: transmission protocol to be ued (`UDP` or `TCP`)
    If None take value from `HARVEST_PROTOCOL` environemental variable, or use default `UDP` if not set.

    :return: tuple with used (host, port, message)
    '''
    import os,socket,copy,random,time

    version=3

    if host is None:
        if 'HARVEST_HOST' in os.environ:
            host=os.environ['HARVEST_HOST']
        else:
            host='gadb-harvest.ddns.net'

    if protocol is None:
        if 'HARVEST_PROTOCOL' in os.environ:
            protocol=os.environ['HARVEST_PROTOCOL']
        else:
            protocol='UDP'

    if port is None:
        if 'HARVEST_PORT' in os.environ:
            port=int(os.environ['HARVEST_PORT'])
        elif protocol=='UDP':
            port=32000
        else:
            port=31000

    if verbose is None:
        if 'HARVEST_VERBOSE' in os.environ:
            verbose=int(os.environ['HARVEST_VERBOSE'])
        else:
            verbose=0

    payload=copy.deepcopy(payload)
    payload['_user']=os.environ['USER']
    payload['_hostname']=socket.gethostname()
    payload['_workdir']=os.getcwd()

    if '_tag' not in payload:
        if tag is None:
            if 'HARVEST_TAG' in os.environ:
                tag=os.environ['HARVEST_TAG']
            else:
                tag=''
        payload['_tag']=tag

    message = "%d:%s:%s"%(version,table,_data_2_message(payload))

    #UDP connection with application level fragmentation
    MTU=1450
    if protocol in ['UDP']:
        if len(message)<MTU:
            try:
                if protocol=='UDP':
                    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
                    sock.sendto(message, (host,port))
                else:
                    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                    sock.connect((host,port))
                    sock.sendall(message)
                    sock.close()
                if verbose:
                    print("%s:%d --%s--[%3.3f]-> %s"%(host,port,protocol,len(message)*1./MTU,message))
            except Exception as _excp:
                if verbose:
                    print(repr(_excp))

        else:
            fmt="&%06d&%03d&%03d&"
            n=MTU-len(fmt%(0,0,0))
            split_message=[message[x:x+n] for x in range(0,len(message),n)]
            ID=(random.randint(0,10**len(str(id(n))))+id(n))/999999
            for k,message in enumerate(split_message):
                message = (fmt+'%s')%(ID,k,len(split_message),message)
                try:
                    if protocol=='UDP':
                        sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
                        sock.sendto(message, (host,port))
                    else:
                        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                        sock.connect((host,port))
                        sock.sendall(message)
                        sock.close()
                    time.sleep(0.01)
                    if verbose:
                        print("%s:%d --%s--[%3.3f]-> %s"%(host,port,protocol,len(message)*1./MTU,message))
                except Exception as _excp:
                    if verbose:
                        print(repr(_excp))

    #TCP connection with fragmentation at TCP layer
    else:
        try:
            sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            sock.connect((host,port))
#            sock.setsockopt(socket.SOL_SOCKET, socket.SO_SNDBUF, MTU*256)
            if False:
                sock.sendall(message)
                if verbose:
                    print("%s:%d --%s--> %s"%(host,port,protocol,message))
            else:
                messages=[message[i:i+MTU] for i in range(0, len(message), MTU)]
                while len(messages):
                    #time.sleep(0.01)
                    tmp=messages.pop(0)
                    sock.send(tmp)
                    if verbose:
                        print("%s:%d --%s--> %s"%(host,port,protocol,tmp))
        except Exception as _excp:
            if verbose:
                print(repr(_excp))
        finally:
            try:
                sock.close()
            except:
                pass

    return (host,port,message)

def harvest_nc(filename, entries=None, verbose=False):
    '''
    Function that returns data contained in a NetCDF3 file to be sent by harvest_send

    :param filename: NetCDF3 file

    :param entries: subset of variables to loof for in the NetCDF file

    :param verbose: print payload

    :return: payload dictionary
    '''
    import os,netCDF4

    payload={}
    payload['_harvest_filename']=os.path.abspath(filename)

    nc = netCDF4.Dataset(filename,'r',format='NETCDF3_CLASSIC')
    if entries is None:
        entries=nc.variables.keys()
    for entry in entries:
        if entry in nc.variables.keys():
            try:
                value=nc.variables[entry].getValue()[0]
            except Exception:
                value=nc.variables[entry][:]
            if len(value)==1:
                payload[entry]=value[0]
            else:
                payload[entry]=value
            if verbose:
                print(str(entry), value[0])
    nc.close()

    return payload

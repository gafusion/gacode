# Map numpy to python string
def gapystr_get(s):
   try:
      n = len(s)
   except:
      return str(s,'utf-8')

   u = []
   for i in range(len(s)):
      u.append(str(s[i],'utf-8').strip())
   return u

# Map python to numpy string
def gapystr_set(s, l=10, n=200):
   return [item.ljust(l) for item in s + [''] * n][:n]

try:
   from gacode_ext import *

   __all__ = ['gapystr_get',
              'gapystr_set',
              'tgyro',
              'cgyro',
              'gyro',
              'neo',
              'profiles_gen']

   # Pack extensions into pygacode
   tmp={}
   exec('from gacode_ext import *', tmp)
   __all__.extend(list(tmp.keys()))
except:
   # Not using pygacode
   pass   

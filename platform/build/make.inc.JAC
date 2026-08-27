#---------------------------------------------------
# JAC {Portland Group}  
#---------------------------------------------------

IDENTITY="JET analysis cluster PGI+mpich (outdated)"
# Compilers and flags

FC	= mpif90 -module $(GACODE_ROOT)/modules
F77	= mpif90

FOMP    = -mp
FMATH	= -r8 
FOPT	= -fastsse 
FDEBUG	= -g -Ktrap=fp -Mbounds -Minform=warn

# System math libraries

#LMATH = /usr/lib/liblapack.a /usr/lib/libblas.a
#LMATH = /usr/local/depot/PGI/pgi-7.0/linux86/7.0-7/lib/liblapack.a \
#        /usr/local/depot/PGI/pgi-7.0/linux86/7.0-7/lib/libblas.a
LMATH = /usr/local/depot/PGI/pgi-12.4/linux86/12.4/lib/liblapack.a \
        /usr/local/depot/PGI/pgi-12.4/linux86/12.4/lib/libblas.a

# netCDF 

NETCDF = /usr/local/depot/netcdf-4.0.1/lib/libnetcdf.a
NETCDF_INC = /usr/local/depot/netcdf-4.0.1/include

# Archive 

ARCH = ar cr

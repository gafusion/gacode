#!/bin/sh

# Get rid of the old file
rm -f make_depends make.sources

# Loop over source file directories
SRC=`ls *.[Ff]{,90} */*.[Ff]{,90}`
SRCS=`echo "$SRC" | grep -v shared_modules/events.f90 | grep -vw cray341.f`
#echo "$SRCS"
#exit
echo 'SOURCES = \' >> make.sources
for file in $1 #iterdb_v1.29.f $SRCS # tglf_modules.o iterdb_v1.29.f 
do
  basen=`python -c "import os;print '.'.join(os.path.split('$file')[1].split('.')[:-1])"`
  f=`echo "$SRCS" | grep -w $basen | head -1`
  echo "$f"'\' >> make.sources
  echo "#$f" >> make_depends
  #echo "$patt"
  # Loop over all source files
  #for f in `ls  $patt` #
  #do

    # Get a list of all USE statements not in commented lines
    # Parse out the module name
    # Add .mod to the module name
    # Convert the module names to lower case
    # Sort the modules
    # Pick out the unique module dependencies
    mod_deps=`grep '^[^c!]' $f | cut -d, -f1 | cut -d! -f1 | awk '{print $1 " " $2}' | grep -i "^USE\b" | grep -v "\bnetcdf\b" | grep -v "\bmpi[,]?\b" | awk '{print $2 ".mod"}' | tr '[:upper:]' '[:lower:]' | tr -d "\r" | sort | uniq | grep -vF "mpi.mod" | grep -vF "=.mod" | grep -vF mpi12.mod | grep -vF typesizes.mod | grep -vF p_nfreya_interface.mod | grep -vF xsct.mod | grep -vF zonal_data.mod | grep -vF nfreya_version.mod`
    #echo $mod_deps
    # Get a list of all MODULE statements not in commented lines or are procedures
    # Parse out the module name
    # Add .mod to the module name
    # Convert the module name to lower case
    # Sort the modules
    mods=`grep '^[^c!]' $f | cut -d! -f1 | grep -i "\bMODULE\b" | awk '{print $1 " " $2}' | grep -i '^module\b' | awk '{print $2 ".mod"}' | tr '[:upper:]' '[:lower:]' | tr -d "\r" | sort | uniq | grep -vF procedure.mod`
    #echo "$mods"
    # Get only the external modules
    ext_mods=`(echo "$mods" & echo "$mods" & echo "$mod_deps") | sort | uniq -u`
    #Only make a dependency list if there were any dependencies
    #if [ `echo "$ext_mods" | wc -c` -gt 1 ] 
    #then
      # Get the filename sans extension
      
      
      # Add the dependency to make_depends
      echo "$basen.o:" $ext_mods >> make_depends
      nmod=`echo "$mods" | wc -c`
      if [ $nmod -eq 1 ] 
      then
        printf '\t$(COMPILE) -I$(NETCDF_INC) -DONETWO '"$f"'\n' >> make_depends
      else #if [ $nmod -eq 1 ] 
        for m in $mods
        do
          printf '\t$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo '"$m"' | tr A-Z a-z | cut -d. -f1` '"$f"'\n' >> make_depends
        done
      fi
      # Loop over all modules defined in this file
      for m in $mods
      do

        #Add a dependency line for these modules
        echo "$m:" $ext_mods >> make_depends
        printf '\t$(COMPILE) -I$(NETCDF_INC) -DONETWO '"$f"'\n' >> make_depends
        #exit
      done # Modules
    
    #fi
  #done # Files
  #exit
done # Directories

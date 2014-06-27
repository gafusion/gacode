#!/bin/sh
echo Processing $1 to determine dependencies
f=$1 # Source file 
outfile=$2 # Place to output the dependencies
basen=`python -s -S -E -c "import os;print '.'.join(os.path.split('$f')[1].split('.')[:-1])"`
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

# Add the dependency to the output make include file
echo "$basen.o:" $ext_mods "$f" > $outfile
nmod=`echo "$mods" | wc -c`
if [ $nmod -eq 1 ] 
then
  printf '\t$(COMPILE) -I$(NETCDF_INC) -DONETWO '"$f"'\n' >> $outfile
else 
  for m in $mods
  do
    printf '\t$(COMP_MOD) "$(COMPILE) -I$(NETCDF_INC) -DONETWO" `echo '"$m"' | tr A-Z a-z | cut -d. -f1` '"$f"'\n' >> $outfile
  done
fi
# Loop over all modules defined in this file
for m in $mods
do

  #Add a dependency line for these modules
  echo "$m:" $ext_mods >> $outfile
  printf '\t$(COMPILE) -I$(NETCDF_INC) -DONETWO '"$f"'\n' >> $outfile
done # Modules

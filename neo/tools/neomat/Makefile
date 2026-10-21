# Accuracy tests of the NEO energy-basis matrices.
# "old" is neo_energy_grid.f90 before the quadrature rewrite.
OLD = 1880d196
SRC = ../../src
LIB = $(SRC)/neo_globals.f90
DEP = $(SRC)/neo_error.f90 $(SRC)/neo_compute_fcoll.f90
R8  = gfortran -O2 -fdefault-real-8 -fdefault-double-8
R16 = gfortran -O2 -freal-4-real-16

all: basis coll

old_energy_grid.f90:
	git show $(OLD):neo/src/neo_energy_grid.f90 > $@

# streaming/source arrays: python3 check_basis.py
basis: old_energy_grid.f90
	mkdir -p b8 b16 o8 o16
	$(R8)  -Jb8  $(LIB) $(SRC)/neo_energy_grid.f90 $(DEP) test_basis.f90 -o new8
	$(R16) -Jb16 $(LIB) $(SRC)/neo_energy_grid.f90 $(DEP) test_basis.f90 -o new16
	$(R8)  -Jo8  $(LIB) old_energy_grid.f90 $(DEP) test_basis.f90 -o old8
	$(R16) -Jo16 $(LIB) old_energy_grid.f90 $(DEP) test_basis.f90 -o old16

# collision matrices: python3 check_coll.py
coll: old_energy_grid.f90
	mkdir -p b8 b16 o8 o16
	$(R8)  -Jb8  $(LIB) $(SRC)/neo_energy_grid.f90 $(DEP) test_coll.f90 -o cnew8
	$(R16) -Jb16 $(LIB) $(SRC)/neo_energy_grid.f90 $(DEP) test_coll.f90 -o cnew16
	$(R8)  -Jo8  $(LIB) old_energy_grid.f90 $(DEP) test_coll.f90 -o cold8
	$(R16) -Jo16 $(LIB) old_energy_grid.f90 $(DEP) test_coll.f90 -o cold16

# old matrix formulas and fcoll vs mpmath: check_matrix.py, check_fcoll.py
drv:
	$(R8)  drv.f90 $(SRC)/neo_compute_fcoll.f90 -o drv8
	$(R16) drv.f90 $(SRC)/neo_compute_fcoll.f90 -o drv16
	$(R8)  fcoll_dump.f90 $(SRC)/neo_compute_fcoll.f90 -o fcoll8
	$(R16) fcoll_dump.f90 $(SRC)/neo_compute_fcoll.f90 -o fcoll16

clean:
	rm -rf b8 b16 o8 o16 new8 new16 old8 old16 cnew8 cnew16 cold8 cold16 \
	       drv8 drv16 fcoll8 fcoll16 old_energy_grid.f90 fcache_*.pkl \
	       out_*.txt __pycache__ *.mod *~

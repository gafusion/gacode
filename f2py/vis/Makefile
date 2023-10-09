all:
	python -m numpy.f2py -c -m vis vis.f90 --f90flags="-fopenmp" -lgomp

clean:
	rm -f *so

integrate: customPolarSolverMcml.pyx
	rm __init__.py
	python setup.py build_ext --inplace
	echo "from customPolarSolverMcml import polarPoissonIntegralSolver" > "__init__.py"

clean:
	rm -r __init__.py* *~ *.c *.so build/ 
	touch __init__.py

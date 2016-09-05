integrate: customCartesianSolverMcxyz.pyx
	rm __init__.py
	python setup.py build_ext --inplace
	echo "from customCartesianSolverMcxyz import cartPoissonIntegralSolver" > "__init__.py"

clean:
	rm -r __init__.py* *~ *.c *.so build/ ; touch __init__.py

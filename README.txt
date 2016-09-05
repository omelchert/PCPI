1. INTRODUCTION

PyPCPI -- A python module for optoacoustic signal generation via [P]olar 
          [C]onvolution and [P]oisson [I]ntegral evaluation

PyPCPI describes a software tool for the numerical calculation of
optoacoustic signals for layered homogeneous media (LHoM) and optically
heterogeneous media (OHeM). In case of LHoM it is designed to postprocess the
material response to an infinitely narrow photon beam obtained via the popular
MCML C code (see http://omlc.org/software/mc/). Therefore, working in
polar coordinates, it implements the functionality to obtain the response to a
spatially extended irradiation source profile via polar convolution and allows
to compute the optoacoustic signal observed by a pointlike detector by solving
the optoacoustic Poission integral.  In case of OHeM it postprocesses the
material response to a spatially extended photon beam obtained via the
well-established MCXYZ C code (see http://omlc.org/software/mc/),
providing an OA Poisson integral solver in cartesian coordinates.  PyPCPI
is implemented using Python following a functional programming paradigm.


2. DEPENDENCIES

The software was developed and tested under OS X Yosemite (Version: 10.10.3)
but should be able to run on any system that features the necessary version 
of Python and has all dependency modules available.

Python -- Version 2.7.6 or higher
numpy  -- Version 1.8.0rc1 or higher
scipy  -- Version 0.13.0b1 or higher 
Cython -- Version 0.15.1 or higher


3. CONTENT

README.txt              -- this readme document
LICENSE.txt             -- BSD 3-Clause License file

PyPCPI/                 -- PyPCPI software module
     __init__.py
     layeredMedia
         __init__.py
         dataIO
             __init__.py
             gpl.py
             mat.py
             mcmlio.py
             misc.py
             npz.py
         poissonIntegralSolver
             __init__.py
             acousticObservables.py
             poissonIntegral_cython
                 __init__.py
                 customPolarSolverMcml.pyx
                 Makefile
                 setup.py
             test
                 test_acousticObservables.py
         polarConvolution
             __init__.py
             convolveRadiallySymmetricFunctions.py
             hankelTransform.py
             irradiationSourceProfile.py
             test
                 test_convolveRadiallySymmetricFunctions.py
                 test_hankelTransform.py
                 test_irradiationSourceProfile.py
     signalPostProcessing
         __init__.py
         acousticAttenuation.py
         finitePulse.py
     voxelizedMedia
         __init__.py
         dataIO
             __init__.py
             gpl.py
             mat.py
             mcxyzio.py
             npz.py
         poissonIntegralSolver
             __init__.py
             acousticObservables.py
             poissonIntegral_cython
                 __init__.py
                 customCartesianSolverMcxyz.pyx
                 Makefile
                 setup.py
             test
                 test_acousticObservables.py
         pureAbsorber
             __init__.py
             irradiationSourceProfile.py
             sourceVolume.py
             test
                 test_beamPropagation.py
                 test_pureAbs.py

exampleUsage/                   -- suit of functional tests for PyPCPI
     mcmlExtension
         inputData
             2Layer_mcml.mci
             2Layer_mcml.mco
         mcml_useCase1_basicExamples
             dataUseCase11
                 Wrz_2Layer_mcml_ISPGauss_R00.15.npz
             dataUseCase12
                 GP
                     FIGS
                         useCase12_Wrz.eps
                         useCase12_Wrz.pdf
                     useCase12_Wrz.gpl
                 useCase12_Wrz.prof
             dataUseCase13
                 GP
                     FIGS
                         useCase13_oaSignal.eps
                     useCase13_oaSignal.gpl
                 useCase13_oaSig_rD0.0_zD-2.0.dat
             runUseCases.sh
             useCase11_main_convolve.py
             useCase12_main_yieldGnuplotOutput.py
             useCase13_main_oaSIgnal.py
         mcml_useCase2_illustrateISPs
             dataUseCase21
                 GP
                     FIGS
                         illustrateISPs_Wrz.eps
                     illustrateISPs_Wrz.gpi
                 Wrz_Donut.npz
                 Wrz_Donut.prof
                 Wrz_flatTop.npz
                 Wrz_flatTop.prof
                 Wrz_Gaussian.npz
                 Wrz_Gaussian.prof
             dataUseCase22
                 GP
                     FIGS
                         oaSignal_differentISPs.eps
                     oaSignal_differentISPs.gpi
                 oaSignal_Donut.dat
                 oaSignal_flatTop.dat
                 oaSignal_Gaussian.dat
             dataUseCase23
                 GP
                     FIGS
                         mcml_oaSig_donut.eps
                     oaSig_Donut.gpi
                 oaSignal_Donut_zD-0.5.dat
                 oaSignal_Donut_zD-1.0.dat
                 oaSignal_Donut_zD-2.0.dat
                 oaSignal_Donut_zD-4.0.dat
                 oaSignal_Donut_zD-8.0.dat
             runUseCase.sh
             useCase21_main_illustrateISPs.py
             useCase22_main_axialOASignal.py
             useCase23_main_finitePulse.py
     mcxyzExtension
         inputData
             skinvessel_F.bin
             skinvessel_H.mci
             skinvessel_T.bin
         mcxyz_useCase1_skinvessel
             dataUseCase11
                 GP
                     FIGS
                         mcxyz_skinvessel_Wxyz.eps
                     skinvessel_Wxyz.gpi
                     skinvessel_Wxyz2.gpi
                 Wxyz_xScan.dat
                 Wxyz_yScan.dat
             dataUseCase12
                 GP
                     FIGS
                         mcxyz_oaSig_zScan_yScan.eps
                     oaSig_zScan_yScan.gpi
                 oaSig_yD0.01.dat
                 oaSig_yD0.015.dat
                 oaSig_yD0.02.dat
                 oaSig_yD0.025.dat
                 oaSig_yD0.03.dat
                 oaSig_yD0.035.dat
                 oaSig_yD0.04.dat
                 oaSig_yD0.045.dat
                 oaSig_yD0.05.dat
                 oaSig_zD-0.05.dat
                 oaSig_zD-0.1.dat
                 oaSig_zD-0.2.dat
                 oaSig_zD-0.4.dat
                 oaSig_zD-0.8.dat
                 oaSig_zD-1.6.dat
                 oaSig_zD-3.2.dat
             runUseCases.sh
             useCase11_main_WxyzSlice.py
             useCase12_main_OASignalOffAxis.py
     pureAbsorber
         pureAbs_useCase1_PSG1996
             dataUseCase1
                 GP
                     FIGS
                         pureAbs_OASignal.eps
                         pureAbs_Wxyz.eps
                     OASignal.gpi
                     Wxyz.gpi
                 oaSignal_fta0.025.dat
                 oaSignal_fta0.05.dat
                 oaSignal_fta0.075.dat
                 oaSignal_fta0.1.dat
                 oaSignal_fta0.194.dat
                 Wxy0z.prof
             pureAbs_useCase1_main.py
             runUseCase.sh
         pureAbs_useCase2_PS2000
             dataUseCase21
                 GP
                     FIGS
                         oaSignals_fig8.eps
                     OASignal_fig8.gpi
                 L532nm_ma80.0_ftd0.06.dat
                 L532nm_ma80.0_ftd0.3.dat
                 L550nm_ma13.0_ftd0.06.dat
                 L550nm_ma13.0_ftd0.3.dat
             dataUseCase22
                 GP
                     FIGS
                         oaSignals.eps
                     OASignal_fig6.gpi
                 oaSignal_doubleLayer.dat
             pureAbs_useCase21_fig8_main.py
             pureAbs_useCase22_fig6_main.py
             runUseCase.sh


4. LICENSE

BSD 3-Clause License


5. ACKNOWLEDGEMENTS

O. Melchert acknowledges support from the VolkswagenStiftung within the
Nieders\"achsisches Vorab program in the framework of the project Hybrid
Numerical Optics. 


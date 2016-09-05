""" FILE: acousticAttenuation.py

Module implementing functions to obtain an optoacoustic signal that features
characteristic changes due to acoustic attenuation in a lossy medium.
See, e.g. 

"Coupling 3D Monte Carlo light transport in optically heterogeneous tissues
to photoacoustic signal generation",
Jacques, S. L.
Photoacoustics 2 (2014) 137-142

Author: O. Melchert
Date:   25.06.2016 
"""
import numpy as np

def scaleFreqs(pSig,dt,alpha0=0.5,eta=0.2):
        """signal change due to acoustic attenuation

        \param pSig original optoacoustic signal 
        \param dt temporal increment
        \param alpha0 frequency scaling amplitude
        \param eta frequency scaling exponent
        """
        ftSig     = np.fft.rfft(pSig)
        freq      = np.fft.rfftfreq(pSig.shape[0],d=dt)
        ftSig[:] *= alpha0*(freq[:])**eta
        return np.fft.irfft(ftSig)

# EOF: acousticAttenuation.py

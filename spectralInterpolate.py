#!/usr/bin/env python
# 
# claw, 13mar13
# written by Earl Lawrence by idea from Scott van der Wiel

import numpy as np
#import pdb

def spectralInterpolate(Y, t0, maxTurns=10, turnIncrement=0.125, weightNorm=2):
	
	if maxTurns==0:
		return(Y.mean())

	# 1:  DFT
	# Number of visibilities
	nr = Y.size
	# Total length for fft
	nt = nr/turnIncrement
	# Half the number of Fourier frequencies considered
	nf = min(nt/2, maxTurns/turnIncrement)

	if nf==nt/2:
		freqID = range(int(nt))
	else:
		freqID = range(int(nf+1));
		freqID.extend(range(int(nt-nf),int(nt)))

	# Fourier frequencies
	omega = np.array(freqID)/nt

	# Some stuff that creates F
	F = np.fft.fft(Y, n=int(nt))
	F = F[freqID]
	F = np.sqrt(2)*F/nr

	# 2:  Calculate weights
	W = np.abs(F)
	W = W/max(W)
	W = W**weightNorm
	W = W/sum(W)

	# 3:  Rotate FFT per frequency, weight, and sum across frequencies
	rotation = np.exp(2*np.pi*omega*t0*1j)
	Yhat = sum(rotation*(W*F))

	#pdb.set_trace()
	return(Yhat)


def run(maxturns=1):
	testdata = np.random.normal(0, 1, (100,2))
	V = np.zeros((testdata.shape[0]), dtype=complex)
	for i in range(V.size):
		V[i] = testdata[i,0] + 1j*testdata[i,1]
	
	t0 = V.size/2
	return spectralInterpolate(V, t0, maxTurns=maxturns)

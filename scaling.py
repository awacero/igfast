#!/usr/bin/python
#This code outputs the predicted magnitude using PGD and Pd from seismogeodetic data.
import math
import numpy

def PGD(d, r, repi):
        l1=len(d)
        A = -6.687
        B = 1.500
        C = -0.214

        Weight = numpy.exp(-numpy.power(repi,2)/2/numpy.power(min(repi),2))

        W = numpy.zeros([l1,l1])
        for i in range (0, l1):
                W[i,i] = Weight[i]

        G = B+C*(numpy.log10(r))
        b = numpy.log10(d)-A
        M = numpy.linalg.lstsq(numpy.dot(W,G),numpy.dot(W,b),rcond=None)[0]
        UP = numpy.dot(G,M)
        SYN = numpy.power(10,UP+A)
        DAT = d

        IQR = numpy.subtract(*numpy.percentile(DAT-SYN, [75, 25]))
        VR = 100*(1-numpy.sum(numpy.sqrt((DAT-SYN)**2))/numpy.sum(numpy.sqrt((DAT)**2)))

        return(M,VR,IQR)








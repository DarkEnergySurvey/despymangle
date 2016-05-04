#!/usr/bin/env python

import sys
import numpy as np
import pymangle as pym
import pylab

import pyfits as pfy
from astropy import wcs

def make_comp(coadd_fullname, fn_mg, fn_star, fn_bleed, plot_fullname,
              limitx=np.arange(0, 10000, 10),
              limity=np.arange(0, 10000, 10)):
    """ Make comparison between coadd file and mangle outputs """
    mg=pym.Mangle(fn_mg)
    star=pym.Mangle(fn_star)
    bleed=pym.Mangle(fn_bleed)

    h=pfy.open(coadd_fullname)
    w2 = wcs.WCS(h[2].header)

    print "limitx=", limitx
    print "limity=", limity
    A = np.array([limitx, limity])
    M = np.meshgrid(A[1], A[0])
    print "A[1] =", A[1]
    print "A[0] =", A[0]


    fig = pylab.figure(figsize=(30, 5))

    pylab.subplot(141)
    pylab.pcolormesh(A[0], A[1], h[2].data[M].T, 
               vmin=np.min(h[2].data[M]), 
               vmax=np.max(h[2].data[M]))
    pylab.colorbar()
    pylab.title('Coadd weight')

    N = M[0].flatten()
    P = M[1].flatten()
    print "N = ", N
    print "P = ", P
    PN = np.array([P, N]).T
    print "PN =", PN
    print type(PN)
    print PN.shape

    B = w2.wcs_pix2world(np.array([P, N]).T, 1)
    print B

    val = mg.weight(B[:,0], B[:,1])
    val2 = star.weight(B[:,0], B[:,1])
    val3 = bleed.weight(B[:,0], B[:,1])
    val = val*(1-val2)*(1-val3)
    D = val.reshape((len(A[1]), len(A[0])), order='F')

    pylab.subplot(142)
    pylab.pcolormesh(A[0], A[1], D, 
               vmin=np.min(h[2].data[M]), 
               vmax=np.max(h[2].data[M]))
    pylab.colorbar()
    pylab.title('Mangle weight')


    pylab.subplot(143)
    pylab.pcolormesh(A[0], A[1], (h[2].data[M].T-D), 
               vmin=-np.max(h[2].data[M])/5., 
               vmax=np.max(h[2].data[M])/5.)
    pylab.colorbar()
    pylab.title('Difference (Coadd - Mangle)')

    ng = np.where(D==0)
    rap = h[2].data[M].T/D
    rap[ng] = -1
    rap = np.ma.masked_values(rap, -1)

    pylab.subplot(144)
    pylab.pcolormesh(A[0], A[1], rap, vmin=0.9, vmax=1.1)
    pylab.colorbar()
    pylab.title('Ratio (Coadd / Mangle)')


    pylab.savefig(plot_fullname)

    h.close

    del(mg, D, val, A, N ,P, rap)


import numpy as np
import matplotlib
matplotlib.use('Agg')
import pylab
import fitsio
import pymangle as pym

from despyastro import wcsutil
import despymangle.EmptyMangle as em

def make_comp(coadd_fullname, fn_mg, fn_star, fn_bleed, plot_fullname,
              limitx=np.arange(0, 10000, 10),
              limity=np.arange(0, 10000, 10)):
    """ Make comparison between coadd file and mangle outputs """
    f = open(fn_mg, 'r')
    line = f.readline()
    if line.startswith('0 polygons'):
        mg = em.Empty_Mangle()
    else:
        mg = pym.Mangle(fn_mg)

    f = open(fn_star, 'r')
    line = f.readline()
    if line.startswith('0 polygons'):
        star = em.Empty_Mangle()
    else:
        star = pym.Mangle(fn_star)

    f = open(fn_bleed, 'r')
    line = f.readline()
    if line.startswith('0 polygons'):
        bleed = em.Empty_Mangle()
    else:
        bleed = pym.Mangle(fn_bleed)

    data2, h2 = fitsio.read(coadd_fullname, ext='WGT', header=True)
    w2 = wcsutil.WCS(h2)

    A = np.array([limitx, limity])
    M = np.meshgrid(A[1], A[0])

    _ = pylab.figure(figsize=(30, 5))

    pylab.subplot(141)
    pylab.axis([np.min(limitx), np.max(limitx), np.min(limity), np.max(limity)])
    pylab.pcolormesh(A[0], A[1], data2[M].T,
                     vmin=np.min(data2[M]),
                     vmax=np.max(data2[M]))
    pylab.colorbar()
    pylab.title('Coadd weight')

    N = M[0].flatten()
    P = M[1].flatten()

    ra, dec = w2.image2sky(P, N)

    val = mg.weight(ra, dec)
    val2 = star.weight(ra, dec)
    val3 = bleed.weight(ra, dec)
    val = val * (1 - val2) * (1 - val3)
    D = val.reshape((len(A[1]), len(A[0])), order='F')

    pylab.subplot(142)
    pylab.axis([np.min(limitx), np.max(limitx), np.min(limity), np.max(limity)])
    pylab.pcolormesh(A[0], A[1], D,
                     vmin=np.min(data2[M]),
                     vmax=np.max(data2[M]))
    pylab.colorbar()
    pylab.title('Mangle weight')


    pylab.subplot(143)
    pylab.axis([np.min(limitx), np.max(limitx), np.min(limity), np.max(limity)])
    pylab.pcolormesh(A[0], A[1], (data2[M].T-D),
                     vmin=-np.max(data2[M])/5.,
                     vmax=np.max(data2[M])/5.)
    pylab.colorbar()
    pylab.title('Difference (Coadd - Mangle)')

    ng = np.where(D == 0)
    rap = data2[M].T/D
    rap[ng] = -1
    rap = np.ma.masked_values(rap, -1)

    pylab.subplot(144)
    pylab.axis([np.min(limitx), np.max(limitx), np.min(limity), np.max(limity)])
    pylab.pcolormesh(A[0], A[1], rap, vmin=0.95, vmax=1.05)
    pylab.colorbar()
    pylab.title('Ratio (Coadd / Mangle)')


    pylab.savefig(plot_fullname)

    del(mg, star, bleed, D, val, A, N, P, rap)

import numpy as np
import os
import shutil
import shlex
##import set_basic_mangle_params_ABLTEST as par
import matplotlib.pyplot as plt
import despydb
import sys
from despyastro import wcsutil
import subprocess



def runcmd(cmd, manglebindir=None, log=None):
    """ Run a command redirecting stdout/stderr to log file if requested """

    if manglebindir:
        cmd = manglebindir + '/' + cmd

    print "\n\ncmd =", cmd
    sys.stdout.flush()

    mystdpipe = sys.stdout 
    if log:
        mystdpipe = open(log, 'w', 0)

    try:
        process = subprocess.Popen(shlex.split(cmd),
                                        shell=False,
                                        stdout=mystdpipe,
                                        stderr=subprocess.STDOUT)
        process.wait()
    except:
        (extype, exvalue, _) = sys.exc_info()
        print "********************"
        print "Unexpected error: %s - %s" % (extype, exvalue)
        print "cmd> %s" % cmd
        print "Probably could not find %s in path" % shlex.split(cmd)[0]
        print "Check for mispelled execname in cmd or"
        print "    make sure that the corresponding eups package is in the metapackage "
        print "    and it sets up the path correctly"
        raise

    if log:
        close(mystdpipe)

    if process.returncode != 0:
        raise ValueError("non-zero exit code")

    sys.stdout.flush()
    return process.returncode



#def int2str5(i):
    #if i <10:
    #    res='0000%d'%i
    #elif i<100:
    #    res='000%d'%i
    #elif i<1000:
    #    res='00%d'%i
    #elif i<10000:
    #    res='0%d'%i
    #elif i<100000:
    #    res='%d'%i
#    res = '%05d' % i
#
#    return res






def mean_ra(rr1,rr2):
    if rr1<180:
        rr1+=360
    if rr2<180:
        rr2+=360
    return  np.mod((rr1+rr2)/2.0  ,360  )

def compute_corners_and_middles_ampAB(indict, i, b1, b2, ccdnum):

#    
#                      x-axis#
#
#       (RA4,DEC4)                   (RA3,DEC3)
#       Corner 4 +-----------------+ Corner 3
#     (1,NAXIS2) |+---------------+| (NAXIS1,NAXIS2)
#                ||d      g      c||
#                ||               ||
#                ||               ||
#                ||               ||
#                ||               ||
#                ||               ||
#                ||h             f||
#                ||               ||   y-axis
#                ||               ||
#                ||               ||
#                ||               ||
#                ||               ||
#                ||               ||
#                ||a      e      b||
#      (RA1,DEC1)|+---------------+| (RA2,DEC2)
#       Corner 1 +-----------------+ Corner 2
#          (1,1)                     (NAXIS1,1)
# with
#a=(1+b1, 1+b2)
#b=(n1-b1, 1+b2)
#c=(n1-b1, n2-b2)
#d=(1+b1, n2-b2)
#e=(n1/2, 1+b2)
#f=(n1-b1, n2/2 )
#g=(n1/2, n2-b)
#h=(1+b,n2/2 )

    
    hdr = {}
    hdr['crval1'] = float(indict['CRVAL1'][i])
    hdr['crval2'] = float(indict['CRVAL2'][i])
    hdr['crpix1'] = float(indict['CRPIX1'][i])
    hdr['crpix2'] = float(indict['CRPIX2'][i])
    hdr['ctype1'] = indict['CTYPE1'][i]
    hdr['ctype2'] = indict['CTYPE2'][i]
    
    # Real dictionary, no need to index to [0]
    hdr['naxis1'] = float(indict['NAXIS1'][0])
    hdr['naxis2'] = float(indict['NAXIS2'][0])
    
    hdr['cd1_1'] = float(indict['CD1_1'][i])
    hdr['cd1_2'] = float(indict['CD1_2'][i])
    hdr['cd2_1'] = float(indict['CD2_1'][i])
    hdr['cd2_2'] = float(indict['CD2_2'][i])
    
    hdr['pv1_0'] = float(indict['PV1_0'][i])
    hdr['pv1_1'] = float(indict['PV1_1'][i])
    hdr['pv1_2'] = float(indict['PV1_2'][i])
    hdr['pv1_3'] = 0                          # Felipe asks: Why hardcode 0 and not use indict[pv1_3]?
    hdr['pv1_4'] = float(indict['PV1_4'][i])
    hdr['pv1_5'] = float(indict['PV1_5'][i])
    hdr['pv1_6'] = float(indict['PV1_6'][i])
    hdr['pv1_7'] = float(indict['PV1_7'][i])
    hdr['pv1_8'] = float(indict['PV1_8'][i])
    hdr['pv1_9'] = float(indict['PV1_9'][i])
    hdr['pv1_10'] = float(indict['PV1_10'][i])
    
    hdr['pv2_0'] = float(indict['PV2_0'][i])
    hdr['pv2_1'] = float(indict['PV2_1'][i])
    hdr['pv2_2'] = float(indict['PV2_2'][i])
    hdr['pv2_3'] = 0                           # Felipe asks: Why hardcore 0 and not use indict[pv2_3]?
    hdr['pv2_4'] = float(indict['PV2_4'][i])
    hdr['pv2_5'] = float(indict['PV2_5'][i])
    hdr['pv2_6'] = float(indict['PV2_6'][i])
    hdr['pv2_7'] = float(indict['PV2_7'][i])
    hdr['pv2_8'] = float(indict['PV2_8'][i])
    hdr['pv2_9'] = float(indict['PV2_9'][i])
    hdr['pv2_10'] = float(indict['PV2_10'][i])
    
    
    # Short cut for less typing
    nx = hdr['naxis1']
    ny = hdr['naxis2']
    print 'nx, ny =', nx, ny
 
    wcs = wcsutil.WCS(hdr)
    ra_a,dec_a = wcs.image2sky(1+b1, 1+b2)
    ra_b,dec_b = wcs.image2sky( ( nx/2+ (1+b1))/2,1+b2)
    ra_c,dec_c = wcs.image2sky(nx/2,1+b2)
    ra_d,dec_d = wcs.image2sky(nx/2, (ny-b2 +(1+b2)  )/2    )
    ra_e, dec_e= wcs.image2sky(nx/2, ny-b2 )
    ra_f, dec_f= wcs.image2sky(( nx/2+ (1+b1))/2 , ny-b2 )
    ra_g, dec_g= wcs.image2sky(1+b1 , ny-b2 )
    ra_h, dec_h= wcs.image2sky(1+b1, (ny-b2 +(1+b2)  )/2 )


    ra_aa,dec_aa = wcs.image2sky(1+nx/2, 1+b2)
    ra_bb,dec_bb = wcs.image2sky( ( nx-b1 + (1+nx/2))/2.0,1+b2)
    ra_cc,dec_cc = wcs.image2sky(nx-b1,1+b2)
    ra_dd,dec_dd = wcs.image2sky(nx-b1, (ny-b2 +(1+b2)  )/2    )
    ra_ee, dec_ee= wcs.image2sky(nx-b1, ny-b2 )
    ra_ff, dec_ff= wcs.image2sky((nx-b1 + (1+nx/2))/2 , ny-b2 )
    ra_gg, dec_gg= wcs.image2sky(1+nx/2 , ny-b2 )
    ra_hh, dec_hh= wcs.image2sky(1+nx/2, (ny-b2 +(1+b2)  )/2 )

    edges1bis=c2e(ra_a, dec_a,ra_c, dec_c,ra_e, dec_e,ra_g, dec_g)[0]
    edges2bis=c2e(ra_aa, dec_aa,ra_cc, dec_cc,ra_ee, dec_ee,ra_gg, dec_gg)[0]
    edge1=(ra_a, dec_a, ra_b, dec_b, ra_c, dec_c, ra_d, dec_d, ra_e, dec_e,ra_f, dec_f, ra_g, dec_g, ra_h, dec_h)
    edge2=(ra_aa, dec_aa, ra_bb, dec_bb, ra_cc, dec_cc, ra_dd, dec_dd, ra_ee, dec_ee,ra_ff, dec_ff, ra_gg, dec_gg, ra_hh, dec_hh)

    if ccdnum <32:
        ampA=edges1bis
        ampB=edges2bis
    else:
        ampA=edges2bis
        ampB=edges1bis

    return ampA, ampB










def compute_corners_and_middles(indict, i, b1, b2):

#    
#                      x-axis#
#
#       (RA4,DEC4)                   (RA3,DEC3)
#       Corner 4 +-----------------+ Corner 3
#     (1,NAXIS2) |+---------------+| (NAXIS1,NAXIS2)
#                ||d      g      c||
#                ||               ||
#                ||               ||
#                ||               ||
#                ||               ||
#                ||               ||
#                ||h             f||
#                ||               ||   y-axis
#                ||               ||
#                ||               ||
#                ||               ||
#                ||               ||
#                ||               ||
#                ||a      e      b||
#      (RA1,DEC1)|+---------------+| (RA2,DEC2)
#       Corner 1 +-----------------+ Corner 2
#          (1,1)                     (NAXIS1,1)
# with
#a=(1+b1, 1+b2)
#b=(n1-b1, 1+b2)
#c=(n1-b1, n2-b2)
#d=(1+b1, n2-b2)
#e=(n1/2, 1+b2)
#f=(n1-b1, n2/2 )
#g=(n1/2, n2-b)
#h=(1+b,n2/2 )

    
    hdr = {}
    hdr['crval1'] = float(indict['CRVAL1'][i])
    hdr['crval2'] = float(indict['CRVAL2'][i])
    hdr['crpix1'] = float(indict['CRPIX1'][i])
    hdr['crpix2'] = float(indict['CRPIX2'][i])
    hdr['ctype1'] = indict['CTYPE1'][i]
    hdr['ctype2'] = indict['CTYPE2'][i]
    
    # Real dictionary, no need to index to [0]
    hdr['naxis1'] = float(indict['NAXIS1'][0])
    hdr['naxis2'] = float(indict['NAXIS2'][0])
    
    hdr['cd1_1'] = float(indict['CD1_1'][i])
    hdr['cd1_2'] = float(indict['CD1_2'][i])
    hdr['cd2_1'] = float(indict['CD2_1'][i])
    hdr['cd2_2'] = float(indict['CD2_2'][i])
    
    hdr['pv1_0'] = float(indict['PV1_0'][i])
    hdr['pv1_1'] = float(indict['PV1_1'][i])
    hdr['pv1_2'] = float(indict['PV1_2'][i])
    hdr['pv1_3'] = 0                          # Felipe asks: Why hardcode 0 and not use indict[pv1_3]?
    hdr['pv1_4'] = float(indict['PV1_4'][i])
    hdr['pv1_5'] = float(indict['PV1_5'][i])
    hdr['pv1_6'] = float(indict['PV1_6'][i])
    hdr['pv1_7'] = float(indict['PV1_7'][i])
    hdr['pv1_8'] = float(indict['PV1_8'][i])
    hdr['pv1_9'] = float(indict['PV1_9'][i])
    hdr['pv1_10'] = float(indict['PV1_10'][i])
    
    hdr['pv2_0'] = float(indict['PV2_0'][i])
    hdr['pv2_1'] = float(indict['PV2_1'][i])
    hdr['pv2_2'] = float(indict['PV2_2'][i])
    hdr['pv2_3'] = 0                           # Felipe asks: Why hardcore 0 and not use indict[pv2_3]?
    hdr['pv2_4'] = float(indict['PV2_4'][i])
    hdr['pv2_5'] = float(indict['PV2_5'][i])
    hdr['pv2_6'] = float(indict['PV2_6'][i])
    hdr['pv2_7'] = float(indict['PV2_7'][i])
    hdr['pv2_8'] = float(indict['PV2_8'][i])
    hdr['pv2_9'] = float(indict['PV2_9'][i])
    hdr['pv2_10'] = float(indict['PV2_10'][i])
    
    
    # Short cut for less typing
    nx = hdr['naxis1']
    ny = hdr['naxis2']
    print 'nx, ny =', nx, ny
 
    wcs = wcsutil.WCS(hdr)
    ra_a,dec_a = wcs.image2sky(1+b1, 1+2)
    ra_b,dec_b = wcs.image2sky(nx-b1,1+b2)
    ra_c,dec_c = wcs.image2sky(nx-b1,ny-b2)
    ra_d,dec_d = wcs.image2sky(1+b1, ny-b2)
    ra0 ,dec0  = wcs.image2sky(nx/2.0, ny/2.0)
    ra_e, dec_e= wcs.image2sky(nx/2, 1+b2 )
    ra_f, dec_f= wcs.image2sky(nx-b1, ny/2 )
    ra_g, dec_g= wcs.image2sky(nx/2, ny-b2 )
    ra_h, dec_h= wcs.image2sky(1+b1, ny/2 )


#a=(1+b1, 1+b2)
#b=(n1-b1, 1+b2)
#c=(n1-b1, n2-b2)
#d=(1+b1, n2-b2)
#e=(n1/2, 1+b2)
#f=(n1-b1, n2/2 )
#g=(n1/2, n2-b)
#h=(1+b,n2/2 )

    return (ra_a, dec_a, ra_b, dec_b, ra_c, dec_c, ra_d, dec_d, ra_e, dec_e,ra_f, dec_f, ra_g, dec_g, ra_h, dec_h )


def c2e(ra1, dec1, ra2, dec2, ra3, dec3, ra4,dec4):
    def rameam(rr1,rr2):
        if rr1<180:
            rr1+=360
        if rr2<180:
            rr2+=360
        return  np.mod((rr1+rr2)/2.0  ,360  )

    r1=np.mod(ra1, 360)
    d1=dec1

    r2= rameam(ra1,ra2)
    d2= 0.5*(  dec1+dec2   )

    r3=np.mod(ra2, 360)
    d3=dec2

    r4= rameam(ra2,ra3)
    d4= 0.5*(  dec2+dec3   )

    r5=np.mod(ra3, 360)
    d5=dec3

    r6= rameam(ra3,ra4)
    d6= 0.5*(  dec3+dec4   )

    r7=np.mod(ra4, 360)
    d7=dec4

    r8= rameam(ra4,ra1) 
    d8= 0.5*(  dec4+dec1   )

    s= "%.12f "*16
    a=(r1,d1,r2,d2,r3,d3,r4,d4,r5,d5,r6,d6,r7,d7,r8,d8)

    return a, s%a

def c2v(ra1, dec1, ra2, dec2, ra3, dec3, ra4,dec4):
    s= "%.9f "*8
    a=(ra1,dec1,ra2,dec2,ra3,dec3,ra4,dec4)

    return a, s%a

def create_weighted_poly_from_edges(ra1, dec1, ra2, dec2, ra3, dec3, ra4,dec4, id, w, polyname):
    
    
    g=open('jpol', 'w')
    aa=c2e(ra1, dec1, ra2, dec2, ra3, dec3, ra4,dec4)
    print >>g, "edges %s"%id
    print >>g, aa[1]
    g.close()
    os.system('poly2poly -q  -ie2 jpol jpol1')
  
    os.system(' echo %f >> jw'%w)
    os.system(' weight  -q  -zjw jpol1 %s'%polyname)
    os.system(' rm  jpol jpol1 jw')
    
def create_weighted_poly_from_vertice(ra1, dec1, ra2, dec2, ra3, dec3, ra4,dec4, id, w, polyname):
    
    
    g=open('jpol', 'w')
    aa=c2v(ra1, dec1, ra2, dec2, ra3, dec3, ra4,dec4)
    print >>g, "vertices %s"%id
    print >>g, aa[1]
    g.close()
    os.system('poly2poly -q  -iv4 jpol jpol1')
  
    os.system(' echo %f >> jw'%w)
    os.system(' weight  -q  -zjw jpol1 %s'%polyname)
    os.system(' rm  jpol1 jpol jw')
    
def create_weighted_poly_from_edges16(ra1, dec1, ra2, dec2, ra3, dec3, ra4,dec4,ra5, dec5, ra6, dec6, ra7, dec7, ra8,dec8, id, w, polyname):
    
    
    g=open('jpol', 'w')
    #    aa=c2e(ra1, dec1, ra2, dec2, ra3, dec3, ra4,dec4)
    s= "%.9f "*16

    
    print >>g, "edges %s"%id
    print >>g, s%(ra1, dec1, ra2, dec2, ra3, dec3, ra4,dec4,ra5, dec5, ra6, dec6, ra7, dec7, ra8,dec8)
    g.close()
    os.system('poly2poly  -q  -ie2 jpol jpol1')
    
    os.system(' echo %f >> jw'%w)
    os.system(' weight  -q  -zjw jpol1 %s'%polyname)
    os.system(' rm jpol jpol1 jw')
    


def get_tilename(runid):
    return runid.split('_')[-1]


def get_tileid(tiledir, project, tilename):

    try:
        desdmfile = os.environ["des_services"]
    except KeyError:
        desdmfile = None
    
    section="db-desoper"
    dbh = despydb.desdbi.DesDbi(desdmfile,section)



    # Prepare the query

    query="""select coaddtile_id from coaddtile where  tilename='%s'"""%tilename


    # Get the cursor
    cur=dbh.cursor()
    cur.execute(query)


    item = cur.fetchone()
    return item[0]




#    f=open(tiledir+project+'.tile')
#    tileid=-1
 #   while tileid==-1:
  #      line=f.readline().strip().split(',')
   #     if line[2]==tilename:
    #        tileid=line[0]
    #return tileid

def get_tileid_old(tiledir, project, tilename):

    f=open(tiledir+project+'.tile')
    tileid=-1
    while tileid==-1:
        line=f.readline().strip().split(',')
        if line[2]==tilename:
            tileid=line[0]
    return tileid



def get_neighbor_tolys(project, tileid,tiledir,  outfile, log):

    print "finding tolygons for neighbors..."
    tilefile='%s%s_tiles%s.dpol/%s_tiles%s.%s.pol'%(tiledir, project,par.restag,project,  par.restag, tileid)

    alltilefile='%s%s_tiles%s.pol'%(tiledir, project,par.restag)
    alltolyfile='%s%s_tolys%s.pol'%(tiledir, project,par.restag)

    shutil.copyfile(tilefile, 'jtilefile')
    shutil.copyfile(alltilefile, 'jalltilefile')
    shutil.copyfile(alltolyfile, 'jalltolyfile')

 # make list of id numbers of neighboring tiles by trimming the polygon file containing the tiles with the single tile file for the target tile
    os.system( par.manglebindir+'/rasterize -T -q %s %s %s' %('jtilefile', 'jalltilefile', 'jneigtile'))
    os.system( par.manglebindir+'/poly2poly -q -oi %s %s' %('jneigtile',  'jneightilelist'  ))
    neightilelist=np.unique(np.loadtxt('jneightilelist', skiprows=2))

    os.system('echo %s > jweightedge'%par.edgemaskval)
    os.system('echo 0 > jweight0')

    tolycount=0


    for ntileid in  neightilelist:
        
        ntileid='%d'%ntileid
        tolycount+=1
        tolyfile='%s/%s_tolys%s.dpol/%s_tolys%s.%s.pol'%(tiledir, project, par.restag, project, par.restag, ntileid)


        os.system('cp %s %s'%(tolyfile, 'jtolyfile'))

        #weight neighbors by $edgemaskval, central tile by 0
        if tileid == ntileid:
            os.system( par.manglebindir+'/weight -q -z%s %s %s '%('jweight0',  'jtolyfile',  'jneighborpolys%s'%tolycount  ))
        else:
            os.system( par.manglebindir+'/weight -q -z%s %s %s '%('jweightedge',  'jtolyfile',  'jneighborpolys%s'%tolycount  ))

    
    jneighborpolys=['jneighborpolys%d'%tolycnt for tolycnt in range(1, tolycount+1)]

    jj=''
    for jjj in jneighborpolys:
        jj=jj+jjj+' '
    os.system( par.manglebindir+'/poly2poly -q %s %s >>%s '%(jj, 'jout', log))
    os.system( par.manglebindir+'/snap  '+ par.snap +' %s %s'%('jout', 'jsnap')+ ' >> %s'%log)

    os.system( 'mv %s %s'%( 'jsnap' ,outfile))
    os.system('rm %s'%jj)
    os.system('rm %s %s  %s  %s'%('jtolyfile', 'jneightilelist', 'jneigtile', 'jout'))

   

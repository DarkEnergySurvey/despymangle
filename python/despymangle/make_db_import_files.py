#!/usr/bin/env python

# $Id$
# $Rev::                                  $:  # Revision of last commit.
# $LastChangedBy::                        $:  # Author of last commit.
# $LastChangedDate::                      $:  # Date of last commit.

import despymangle.mangle_db as mdb
import numpy as nm
import numpy as np
import os
import pymangle as pym


def id2indices(Image_tab, ids):
    
    indices = []
    last_digit=[]
    for one_id in ids:
        # drop last digit which means sec A vs B of image 
        img_id = int(str(one_id)[:-1])   
        ldig=int(str(one_id)[-1]) 
        # find the index in the Image_tab where the image with that mangle_img_id exists
        idx = np.where(Image_tab["MANGLE_IMG_ID"] == img_id)[0][0]
        indices.append(idx)
        last_digit.append(ldig)
    return indices, last_digit



def make_csv_files(config,Image_tab,  dbi):

    print config
    if config['outputdir'] is not None:
        outdir = config['outputdir'] + '/'
    fnprefix = outdir + config['molysprefix']
    print '!!!!! ',   config['molysprefix'], outdir




    coadd_object_tab = mdb.get_coadd_object_info(dbi, 
                                                 config['tilename'],
                                                 config['pfw_attempt_id'],
                                                 config['schema'], 
                                                 config['ra_column'],
                                                 config['dec_column'])

    #print coadd_object_tab
    #print coadd_object_tab['RA']
    #print coadd_object_tab['DEC']
    #print coadd_object_tab['ID']


    fndb_pattern = config['outputdir'] + config['dbcsv'] + '_%s.csv'
    fndb_ccdgon = fndb_pattern % ('ccdgon')
    fndb_molygon = fndb_pattern % ('molygon')
    fndb_molygon_ccdgon = fndb_pattern % ('molygon_ccdgon')
    fndb_coadd_object_molygon = fndb_pattern % ('coadd_object_molygon') 

    print 'fndb_ccdgon =', fndb_ccdgon
    print 'fndb_molygon =', fndb_molygon
    print 'fndb_molygon_ccdgon =', fndb_molygon_ccdgon
    print 'fndb_coadd_object_molygon =', fndb_coadd_object_molygon




    #################################
    #### make coadd_object_molygon
    #################################

    #### maybe implement something to check whether there is a mangle mask in this band ?

    ra=coadd_object_tab['RA']
    dec=coadd_object_tab['DEC']

    toly_mask=pym.Mangle(config['poltolys'])
    mangle_mask=pym.Mangle(config['fn_maglims'])
    star_mask=pym.Mangle(config['fn_mask_star'])
    bleed_mask=pym.Mangle(config['fn_mask_bleed'])

    N_obj=len(ra)
    print N_obj




    ##### Check whether object are in the tolys:  since this si in a routine that will be call 5 times, it will be re-checked 5 times whereas once suffice ...

    tol=toly_mask.polyid(ra,dec)
    print 'there are ', len(tol[tol<0]) ,' objects in the tiles which are lost by Mangle'


    import matplotlib.pyplot as plt
    plt.figure(1)
    plt.plot(ra[tol<0], dec[tol<0], 'k.')
    plt.savefig('Object_not_in_tolygon_%s.png'%config['molysprefix'])


    ### Check whether object are  in the geometry mask. 
    A=mangle_mask.polyid_and_weight(ra,dec)   ## A[0] give the id of the polygon (-1 if no polygon), and A[1] give the associated mag_lim
    print A
    print A[1][A[0]<10]

    plt.figure(2)
    plt.plot(ra[A[0]<0], dec[A[0]<0], 'k.')
    plt.savefig('Objects_not_in_geometry_mask_%s.png'%config['molysprefix'])



    ######### Check whether object are in the star_mask

    star=star_mask.polyid(ra,dec)

    plt.figure(3)
    plt.plot(ra[star>=0], dec[star>=0], 'k.')
    plt.savefig('Objects_in_star_mask_%s.png'%config['molysprefix'])


    ######### Check whether object are in the bleed_mask

    bleed=bleed_mask.polyid(ra,dec)

    plt.figure(4)
    plt.plot(ra[bleed>0], dec[bleed>0], 'k.')
    plt.savefig('Objects_in_bleed_mask_%s.png'%config['molysprefix'])



    #### BIT values
    ISNOT_INGEOMETRY=1
    IS_INSTAR=2
    IS_INBLEED=4
    ISNOT_INTOLY=8



    



    basicline='%d,%d,%d'


    g=open(fndb_coadd_object_molygon, 'w')

    for i in range(0,N_obj):
        MOLYGON_ID_X=A[0][i]
        MANGLE_FLAG_X= ISNOT_INGEOMETRY * (A[0][i]<0) + IS_INSTAR *(star[i]>=0) + IS_INBLEED *(  bleed[i]>=0  ) + ISNOT_INTOLY *(tol[i]<0   )
     
        print >>g, basicline%(coadd_object_tab['ID'][i], MOLYGON_ID_X,MANGLE_FLAG_X)
    g.close()










    ##################################
    ##### Make   CCDGON .csv file
    ###################################
    N_images = len(Image_tab['FILENAME'])
    print N_images


    os.system('poly2poly -oi %s toto'%config['fn_ccdgon'] )
    ids=nm.loadtxt('toto', skiprows=2, dtype=nm.int64)


    N_polys=len(ids)
    print N_polys


    RED_IMAGE_FILENAME= Image_tab['FILENAME'][id2indices(Image_tab, ids)[0]]


    os.system('rm toto')

    #print config

    COADD_FILENAME =' something_I_dont_know'
    


    indA=nm.where(id2indices(Image_tab, ids)[1] ==0)
    indB=nm.where(id2indices(Image_tab, ids)[1] ==1)

    os.system('poly2poly -ow %s titi'%config['fn_ccdgon'] )
    weights=nm.loadtxt('titi', skiprows=2)[:,0]

    os.system('rm titi')
    
    basicline='%d,%s,%s,%s,%f,%s'
    g=open(fndb_ccdgon, 'w')

    for i in range(0,N_polys):
        #print i

        if (id2indices(Image_tab, ids)[1][i] ==0):
            AMP='A'
        if (id2indices(Image_tab, ids)[1][i] ==1):
            AMP='B'

        print >>g, basicline%(ids[i], RED_IMAGE_FILENAME[i], COADD_FILENAME,AMP, weights[i],config['pfw_attempt_id'] )
    g.close()





    ###########################
    ### Make MOLYGON .csv file
    ##########################
    

 
   
    os.system('poly2poly -oi  %s_weight.pol   toto'%fnprefix )

    ids=nm.loadtxt('toto', skiprows=2, dtype=nm.int64)
    N_polys=len(ids)
    print N_polys

    COADD_FILENAME =' something_I_dont_know'
    COADDTILE_ID = config['tileid']
    
    ### TOLYGON_ID can probably be removed ??? ABL
    BAND=config['band']
    
    
    NUM_IMAGES=nm.loadtxt('%s/%s.count'%(os.getcwd(),fnprefix), dtype=nm.int32)  ##### ABL: not so sure about the os.getcwd , but loadtxt seem to need absolute paths 
    TOTAL_EXPTIME=nm.loadtxt('%s/%s_EXPTIME.SUM'%(os.getcwd(),fnprefix))
    AREA_STR=nm.loadtxt('%s/%s.area'%(os.getcwd(),fnprefix))
    WAVG_AIRMASS= nm.loadtxt('%s/%s_AIRMASS.WMEAN'%(os.getcwd(),fnprefix))
    WAVG_FWHM=nm.loadtxt('%s/%s_FWHM.WMEAN'%(os.getcwd(),fnprefix))
    MAG_LIMIT=nm.loadtxt('%s/%s.maglims'%(os.getcwd(),fnprefix))
    


    os.system('poly2poly -om %s_weight.pol toto'%fnprefix )

    radec=nm.loadtxt('toto', skiprows=2, usecols=[0,1])

    RA_MID=radec[:,0]
    DEC_MID=radec[:,1]

    basicline='%d,%s,%s,%s,%d,%f,%15.7e,%f,%f,%f,%f,%f,%s'
    i=0
    print ids[i], COADD_FILENAME, COADDTILE_ID, BAND, NUM_IMAGES[i], TOTAL_EXPTIME[i], AREA_STR[i], WAVG_AIRMASS[i], WAVG_FWHM[i], MAG_LIMIT[i], RA_MID[i], DEC_MID[i], config['pfw_attempt_id']

    g=open(fndb_molygon, 'w')

    for i in range(0,N_polys):
      #  print i

        print >>g, basicline%(ids[i], COADD_FILENAME, COADDTILE_ID, BAND, NUM_IMAGES[i], TOTAL_EXPTIME[i], AREA_STR[i], WAVG_AIRMASS[i], WAVG_FWHM[i], MAG_LIMIT[i], RA_MID[i], DEC_MID[i], config['pfw_attempt_id'])
    g.close()




    ###################################
    ### Make molygon_ccdgon .csv file
    ####################################


    #### GET molygons ids


    os.system('poly2poly -oi %s_weight.pol toto'%fnprefix )
    ids=nm.loadtxt('toto', skiprows=2, dtype=nm.int64)
    N_polys=len(ids)
    print N_polys

    
    basicline='%d,%d'
    
    g=open(fndb_molygon_ccdgon, 'w')
    redfile=open('%s/%s.red'%(os.getcwd(),fnprefix))
    #h=open(scripts_dir+'import_in_molygon_ccdgon_table_'+coadd_run+'_'+band+'.dat', 'w')
    for i in range(N_polys):
        aa=redfile.readline().strip().split()
        n=len(aa)-2
        for j in range(0,n):
            print >>g, basicline%(ids[i], nm.int(aa[2+j]))
    g.close()












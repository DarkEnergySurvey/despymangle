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

    mangle_mask=pym.Mangle(config['fn_maglims'])
    star_mask=pym.Mangle(config['fn_mask_star'])
    bleed_mask=pym.Mangle(config['fn_mask_bleed'])


    A=mangle_mask.polyid(ra,dec)
    

    print A
    print A[A<10]






































    ##################################
    ##### Make   CCDGON .csv file
    ###################################
    N_images = len(Image_tab['FILENAME'])
    print N_images


    os.system('poly2poly -oi %s toto'%config['fn_unbalk'] )
    ids=nm.loadtxt('toto', skiprows=2, dtype=nm.int64)


    N_polys=len(ids)
    print N_polys


    RED_IMAGE_FILENAME= Image_tab['FILENAME'][id2indices(Image_tab, ids)[0]]


    os.system('rm toto')

    #print config

    COADD_FILENAME =' something_I_dont_know'
    


    indA=nm.where(id2indices(Image_tab, ids)[1] ==0)
    indB=nm.where(id2indices(Image_tab, ids)[1] ==1)

    os.system('poly2poly -ow %s titi'%config['fn_unbalk'] )
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
    

    if config['outputdir'] is not None:
        outdir = config['outputdir'] + '/'
    fnprefix = outdir + config['molysprefix']
    print '!!!!! ',   config['molysprefix']

    os.system('poly2poly -oi /des001/home/benoitl/y3a1_run/mangle_tests/prodbeta/DES0459-5622/mangle/g/DES0459-5622_r13p68_g_molys_weight.pol toto' )
    ids=nm.loadtxt('toto', skiprows=2, dtype=nm.int64)
    N_polys=len(ids)
    print N_polys

    COADD_FILENAME =' something_I_dont_know'
    COADDTILE_ID = config['tileid']
    
    ### TOLYGON_ID can probably be removed ??? ABL
    BAND=config['band']
    
    NUM_IMAGES=nm.loadtxt('/des001/home/benoitl/y3a1_run/mangle_tests/prodbeta/DES0459-5622/mangle/g/DES0459-5622_r13p68_g_molys.count', dtype=nm.int32)
    TOTAL_EXPTIME=nm.loadtxt('/des001/home/benoitl/y3a1_run/mangle_tests/prodbeta/DES0459-5622/mangle/g/DES0459-5622_r13p68_g_molys_EXPTIME.SUM')
    AREA_STR=nm.loadtxt('/des001/home/benoitl/y3a1_run/mangle_tests/prodbeta/DES0459-5622/mangle/g/DES0459-5622_r13p68_g_molys.area')
    WAVG_AIRMASS= nm.loadtxt('/des001/home/benoitl/y3a1_run/mangle_tests/prodbeta/DES0459-5622/mangle/g/DES0459-5622_r13p68_g_molys_AIRMASS.WMEAN')
    WAVG_FWHM=nm.loadtxt('/des001/home/benoitl/y3a1_run/mangle_tests/prodbeta/DES0459-5622/mangle/g/DES0459-5622_r13p68_g_molys_FWHM.WMEAN')
    MAG_LIMIT=nm.loadtxt('/des001/home/benoitl/y3a1_run/mangle_tests/prodbeta/DES0459-5622/mangle/g/DES0459-5622_r13p68_g_molys.maglims')
    
    os.system('poly2poly -om /des001/home/benoitl/y3a1_run/mangle_tests/prodbeta/DES0459-5622/mangle/g/DES0459-5622_r13p68_g_molys_weight.pol toto' )
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


    os.system('poly2poly -oi /des001/home/benoitl/y3a1_run/mangle_tests/prodbeta/DES0459-5622/mangle/g/DES0459-5622_r13p68_g_molys_weight.pol toto' )
    ids=nm.loadtxt('toto', skiprows=2, dtype=nm.int64)
    N_polys=len(ids)
    print N_polys

    
    basicline='%d,%d'
    
    g=open(fndb_molygon_ccdgon, 'w')
    redfile=open('/des001/home/benoitl/y3a1_run/mangle_tests/prodbeta/DES0459-5622/mangle/g/DES0459-5622_r13p68_g_molys.red')
    #h=open(scripts_dir+'import_in_molygon_ccdgon_table_'+coadd_run+'_'+band+'.dat', 'w')
    for i in range(N_polys):
        aa=redfile.readline().strip().split()
        n=len(aa)-2
        for j in range(0,n):
            print >>g, basicline%(ids[i], nm.int(aa[2+j]))
    g.close()












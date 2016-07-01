#!/usr/bin/env python

# $Id$
# $Rev::                                  $:  # Revision of last commit.
# $LastChangedBy::                        $:  # Author of last commit.
# $LastChangedDate::                      $:  # Date of last commit.

import matplotlib.pyplot as plt
import numpy as np
import os

import pymangle as pym

import despymangle.mangle_db as mdb
import despymangle.mangle_utils as mu
from despymisc import miscutils

######################################################################
def id2indices(Image_tab, ids):
    
    indices  =  []
    last_digit = []
    for one_id in ids:
        # drop last digit which means sec A vs B of image 
        img_id = int(str(one_id)[:-1])   
        ldig = int(str(one_id)[-1]) 
        # find the index in the Image_tab where the image with that mangle_img_id exists
        idx = np.where(Image_tab["MANGLE_IMG_ID"] == img_id)[0][0]
        indices.append(idx)
        last_digit.append(ldig)
    return indices, last_digit


######################################################################
def make_csv_coadd_object_molygon(config, fn_csv, coadd_object_tab):
    """ Gather info and write csv for import into coadd_object_molygon table """

    ra = coadd_object_tab['RA'] 
    dec = coadd_object_tab['DEC']

    toly_mask = pym.Mangle(config['poltolys'])
    mangle_mask = pym.Mangle(config['fn_maglims_pol'])
    star_mask = pym.Mangle(config['fn_mask_star_pol'])
    bleed_mask = pym.Mangle(config['fn_mask_bleed_pol'])

    ##### Check whether object are in the tolys:  since this is in a routine that will be call 5 times, it will be re-checked 5 times whereas once suffice ...
    tol = toly_mask.polyid(ra, dec)
    miscutils.fwdebug_print('There are %d objects in the tiles which are lost by Mangle' % len(tol[tol < 0]))

    ### Check whether objects are in the geometry mask. 
    ## A[0] give the id of the polygon (-1 if no polygon), and A[1] give the associated mag_lim
    A = mangle_mask.polyid_and_weight(ra, dec)
    if miscutils.fwdebug_check(10, "MANGLE_DEBUG"):
        miscutils.fwdebug_print("A = %s" % (A))
        miscutils.fwdebug_print("A[1][A[0] < 10] = %s" % (A[1][A[0] < 10]))

    star = star_mask.polyid(ra, dec)
    bleed = bleed_mask.polyid(ra, dec)

    if config['pltprefix'] is not None:
        plt.figure(1)
        plt.plot(ra[tol < 0], dec[tol < 0], 'k.')
        plt.savefig('%s_Object-not-in-tolygon.png' % config['pltprefix'])

        plt.figure(2)
        plt.plot(ra[A[0] < 0], dec[A[0] < 0], 'k.')
        plt.savefig('%s_Objects_not_in_geometry_masks.png' % config['pltprefix'])

        plt.figure(3)
        plt.plot(ra[star >= 0], dec[star >= 0], 'k.')
        plt.savefig('%s_Objects_in_star_mask.png' % config['pltprefix'])

        plt.figure(4)
        plt.plot(ra[bleed > 0], dec[bleed > 0], 'k.')
        plt.savefig('%s_Objects_in_bleed_mask.png' % config['pltprefix'])


    #### BIT values
    ISNOT_INGEOMETRY = 1
    IS_INSTAR = 2
    IS_INBLEED = 4
    ISNOT_INTOLY = 8
    
    linefmt = '%d,%d,%d'
    g = open(fn_csv, 'w')
    for i in range(0, len(ra)):
        molygon_id_x = A[0][i]
        mangle_flag_x = ISNOT_INGEOMETRY * (A[0][i] < 0) + \
                        IS_INSTAR * (star[i] >= 0) + \
                        IS_INBLEED * (bleed[i] >= 0) + \
                        ISNOT_INTOLY * (tol[i] < 0)
     
        print >>g, linefmt % (coadd_object_tab['ID'][i], molygon_id_x, mangle_flag_x)
    g.close()



######################################################################
def make_csv_ccdgon(config, fn_csv, Image_tab, ra, dec, log):

    jfn_ids = 'jweights_%s_%s' % (config['tileid'], config['band'])
    cmd = 'poly2poly -oi %s %s' % (config['fn_ccdgon_pol'], jfn_ids)
    mu.runcmd(cmd, config['manglebindir'], log)
    if not os.path.exists(jfn_ids):
        raise IOError("Error could not find output file from poly2poly (%s)" % jfn_ids)
    ids = np.loadtxt(jfn_ids, skiprows=2, dtype=np.int64)

    miscutils.fwdebug_print("Number of ccdgon polys = %s" % (len(ids)))

    red_image_filename =  Image_tab['FILENAME'][id2indices(Image_tab, ids)[0]]
    coadd_filename = miscutils.parse_fullname(config['fn_coadd'], miscutils.CU_PARSE_FILENAME)
    
    indA = np.where(id2indices(Image_tab, ids)[1] == 0)
    indB = np.where(id2indices(Image_tab, ids)[1] == 1)

    jfn_weights = 'jweights_%s_%s' % (config['tileid'], config['band'])
    cmd ='poly2poly -ow %s %s' % (config['fn_ccdgon_pol'], jfn_weights)
    mu.runcmd(cmd, config['manglebindir'], log)
    weights = np.loadtxt(jfn_weights, skiprows=2)[:,0]

    id2idx = id2indices(Image_tab, ids)
    
    linefmt = '%d,%s,%s,%s,%f,%s'
    g = open(fn_csv, 'w')
    for i in range(0, len(ids)):
        #print i
        if (id2idx[1][i] == 0):
            amp = 'A'
        elif (id2idx[1][i] == 1):
            amp = 'B'
        else:
            raise KeyError('Invalid amp num for image #%s (%s)' % (i, id2idx[1][i]))
        print >>g, linefmt % (ids[i], red_image_filename[i], coadd_filename, 
                              amp, weights[i], config['pfw_attempt_id'])
    g.close()

    if config['cleanup'] is not None and config['cleanup'].upper() == 'Y':
        os.remove(jfn_weights)


######################################################################
def make_csv_molygon(config, fn_csv, molyids, log):

    jfn_mid = 'jweight-pol-mid_%s_%s' % (config['tileid'], config['band'])
    cmd = 'poly2poly -om %s %s' % (config['fn_molys_weight_pol'], jfn_mid)
    mu.runcmd(cmd, config['manglebindir'], log)
    midradec = np.loadtxt(jfn_mid, skiprows=2, usecols=[0,1])

    ra_mid = midradec[:,0]
    dec_mid = midradec[:,1]
    coadd_filename = miscutils.parse_fullname(config['fn_coadd'], miscutils.CU_PARSE_FILENAME)
    coaddtile_id  =  config['tileid']
    ### TOLYGON_ID can probably be removed ??? ABL
    band = config['band']
    
    num_images = np.loadtxt(config['fn_molys_count'])
    total_exptime = np.loadtxt(config['fn_exptime_sum'])
    area_str = np.loadtxt(config['fn_molys_area'])
    wavg_airmass =  np.loadtxt(config['fn_airmass_wmean'])
    wavg_fwhm = np.loadtxt(config['fn_fwhm_wmean'])
    mag_limit = np.loadtxt(config['fn_molys_maglims'])

    linefmt = '%d,%s,%s,%s,%d,%f,%15.7e,%f,%f,%f,%f,%f,%s'

    if miscutils.fwdebug_check(6, "MANGLE_DEBUG"):
        i = 0
        miscutils.fwdebug_print(linefmt % (molyids[i], coadd_filename, coaddtile_id, 
                                band, num_images[i], total_exptime[i], 
                                area_str[i], wavg_airmass[i], wavg_fwhm[i], 
                                mag_limit[i], ra_mid[i], dec_mid[i], 
                                config['pfw_attempt_id']))

    g = open(fn_csv, 'w')
    for i in range(0, len(molyids)):
        #print i
        print >>g, linefmt % (molyids[i], coadd_filename, coaddtile_id, band, num_images[i], total_exptime[i], area_str[i], wavg_airmass[i], wavg_fwhm[i], mag_limit[i], ra_mid[i], dec_mid[i], config['pfw_attempt_id'])
    g.close()


######################################################################
def make_csv_molygon_ccdgon(config, fn_csv, molyids):

    linefmt = '%d,%d'
    
    g = open(fn_csv, 'w')
    redfh = open(config['fn_molys_red'])
    for i in range(0, len(molyids)):
        aa = redfh.readline().strip().split()
        n = len(aa) - 2
        for j in range(0, n):
            print >>g, linefmt % (molyids[i], np.int(aa[2+j]))
    g.close()

######################################################################
def make_csv_files(config, Image_tab, dbi):
    """ Main function to create csv files for ingestion into db """

    if miscutils.fwdebug_check(10, "MANGLE_DEBUG"):
        miscutils.fwdebug_print("config = %s" % config)

    log = None
    if config['logdir'] is not None:
        log = config['logdir'] + 'make_csv_' + config['tilename'] + '_' + config['band'] +'.log'

    miscutils.fwdebug_print("Number of reduced images: %s" % (len(Image_tab['FILENAME'])))

    ### get information about coadd objects from the DB
    coadd_object_tab = mdb.get_coadd_object_info(dbi, 
                                                 config['tilename'],
                                                 config['pfw_attempt_id'],
                                                 config['schema'], 
                                                 config['ra_column'],
                                                 config['dec_column'])

    
    fndb_pattern  =  config['outputdir'] + config['dbprefix'] + '_%s.csv'

    #### maybe implement something to check whether there is a mangle mask in this band ?

    ra = coadd_object_tab['RA']
    dec = coadd_object_tab['DEC']
    miscutils.fwdebug_print("Number of objects from coadd_object table: %s" % len(ra))

    # get molygon ids
    jfn_ids = 'jweight-pol-ids_%s_%s' % (config['tileid'], config['band'])
    cmd = 'poly2poly -oi %s %s' % (config['fn_molys_weight_pol'], jfn_ids)
    mu.runcmd(cmd, config['manglebindir'], log)
    moly_ids = np.loadtxt(jfn_ids, skiprows=2, dtype=np.int64)
    miscutils.fwdebug_print("Number of ids from moly weight pol = %s" % len(moly_ids))

    if config['cleanup'] is not None and config['cleanup'].upper() == 'Y':
        os.remove(jfn_ids)

    #################################
    #### make coadd_object_molygon .csv file
    #################################
    make_csv_coadd_object_molygon(config, fndb_pattern % 'coadd_object_molygon', coadd_object_tab)

    ##################################
    ##### Make CCDGON .csv file
    ###################################
    make_csv_ccdgon(config, fndb_pattern % 'ccdgon', Image_tab, ra, dec, log)

    ###########################
    ### Make MOLYGON .csv file
    ##########################
    make_csv_molygon(config, fndb_pattern % 'molygon', moly_ids, log)

    ###################################
    ### Make molygon_ccdgon .csv file
    ####################################
    make_csv_molygon_ccdgon(config, fndb_pattern % 'molygon_ccdgon', moly_ids)

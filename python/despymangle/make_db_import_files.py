#!/usr/bin/env python3

# $Id: make_db_import_files.py 45798 2017-06-20 16:01:17Z friedel $
# $Rev:: 45798                            $:  # Revision of last commit.
# $LastChangedBy:: friedel                $:  # Author of last commit.
# $LastChangedDate:: 2017-06-20 11:01:17 #$:  # Date of last commit.


import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import fitsio

import pymangle as pym

import despymangle.mangle_db as mdb
import despymangle.mangle_utils as mu
import despymangle.EmptyMangle as em
from despymisc import miscutils

# used for writing out the systematics
data_types = {'SUM': ('SUM', 'f4'),
              'QSUM': ('QSUM', 'f4'),
              'MIN': ('MIN', 'f4'),
              'MAX': ('MAX', 'f4'),
              'MEDIAN': ('MEDIAN', 'f4'),
              'MEAN': ('MEAN', 'f4'),
              'WMEAN': ('WMEAN', 'f4'),
              'UNCERTAINTY': ('UNCERTAINTY', 'f4')}

######################################################################

def id2skyvar(Image_tab, ids):
    ### return the bood skyvar whether the ccdgons is A or  B

    skyvars = []
    for one_id in ids:
        # drop last digit which means sec A vs B of image
        img_id = int(str(one_id)[:-1])

        is_AB = int(str(one_id)[-1])

        idx = np.where(Image_tab["MANGLE_IMG_ID"] == img_id)[0][0]
        if is_AB == 0:
            skyvars.append(Image_tab["SKYVARA"][idx])
        if is_AB == 1:
            skyvars.append(Image_tab["SKYVARB"][idx])
    return skyvars

######################################################################
def id2indices(Image_tab, ids):

    indices = []
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
def make_syste_masks(keyword, config, Nmolys, Image_tab, molyids, exclude=None, magic_num=0,
                     keys=['SUM', 'MIN', 'MAX', 'MEDIAN', 'MEAN', 'WMEAN']):
    # write out the systematics to a fits file, one HDU for each keyword
    # each of the keys is put into a column, the first column is always the molygonid
    dt = [('MOLYGON_NUMBER', np.int64)]
    for key in keys:
        dt.append(data_types[key])
    tab_exptime = np.zeros(Nmolys, dtype=dt)

    print("\nPROCESSING", keyword)

    f = open(config['fn_molys_red'])
    # for each molygon
    #print config['fn_molys_red']
    for i in range(Nmolys):
        #print i
        line = f.readline().strip()

        ids = line.split()[2:]    ### this is all the ids of the ccdgons in this molygon
        indices, _ = id2indices(Image_tab, ids)  ## indices is the indices in image_Tab which corresponds to the ccdgons
        # get the specific Image_tabs we need.
        if keyword in Image_tab:
            if exclude is not None:
                tmpind = []
                tmpid = []

                for j, val in enumerate(indices):
                    if Image_tab[keyword][j] != exclude:
                        tmpind.append(val)
                        tmpid.append(ids[j])
                indices = tmpind
                ids = tmpid
            temp = Image_tab[keyword][indices]
        else:
            temp = []
        tab_exptime['MOLYGON_NUMBER'][i] = molyids[i]

        if len(temp) != 0:
            # quadtrature sum of the values
            if 'QSUM' in keys and temp is not None:
                tab_exptime['QSUM'][i] = (np.sum(temp ** 2.)) ** 0.5
            # sum of the values
            if 'SUM' in keys and temp is not None:
                tab_exptime['SUM'][i] = np.sum(temp)
                #if keyword == 'FGCM_GRY':
                #    print temp
                #    print " SUM",tab_exptime['SUM'][i]
            # minimum of the values
            if 'MIN' in keys and temp is not None:
                tab_exptime['MIN'][i] = np.min(temp)
            # maximum of the values
            if 'MAX' in keys and temp is not None:
                tab_exptime['MAX'][i] = np.max(temp)
            # median of the values
            if 'MEDIAN' in keys and temp is not None:
                tab_exptime['MEDIAN'][i] = np.median(temp)
            # mean of the values
            if 'MEAN' in keys and temp is not None:
                tab_exptime['MEAN'][i] = np.mean(temp)
            # weighted mean of the values
            if 'WMEAN' in keys and temp is not None:
                skyvars = id2skyvar(Image_tab, ids)
                weights = skyvars * 100 ** ((config['mzpglobal'] - Image_tab['MAG_ZERO'][indices]) / 2.5)
                weights = 1.0 / weights
                tab_exptime['WMEAN'][i] = np.sum(temp * weights) / np.sum(weights)
        # special case of the sky variance
        if 'UNCERTAINTY' in keys:
            skyvars = id2skyvar(Image_tab, ids)
            if skyvars:
                skyvarsp = skyvars * 100 ** ((config['mzpglobal'] - Image_tab['MAG_ZERO'][indices]) / 2.5)
                if 'MIN' in keys:
                    tab_exptime['MIN'][i] = np.min(np.sqrt(skyvarsp))
                if 'MAX' in keys:
                    tab_exptime['MAX'][i] = np.max(np.sqrt(skyvarsp))
                tab_exptime['UNCERTAINTY'][i] = 1. / np.sqrt(np.sum(1. / skyvarsp))

    f.close()

    fits = fitsio.FITS(config['data_file'], 'rw')
    fits.write(tab_exptime, extname=keyword.upper())
    fits.close()



######################################################################
def make_csv_coadd_object_molygon(config, fn_csv, coadd_object_tab):
    """ Gather info and write csv for import into coadd_object_molygon table """

    ra = coadd_object_tab[config['ra_column']]
    dec = coadd_object_tab[config['dec_column']]

    f = open(config['poltolys'], 'r')
    line = f.readline()
    if line.startswith('0 polygons'):
        toly_mask = em.Empty_Mangle()
    else:
        toly_mask = pym.Mangle(config['poltolys'])

    f = open(config['fn_maglims_pol'], 'r')
    line = f.readline()
    if line.startswith('0 polygons'):
        mangle_mask = em.Empty_Mangle()
    else:
        mangle_mask = pym.Mangle(config['fn_maglims_pol'])

    f = open(config['fn_mask_star_pol'], 'r')
    line = f.readline()
    if line.startswith('0 polygons'):
        star_mask = em.Empty_Mangle()
    else:
        star_mask = pym.Mangle(config['fn_mask_star_pol'])

    f = open(config['fn_mask_bleed_pol'], 'r')
    line = f.readline()
    if line.startswith('0 polygons'):
        bleed_mask = em.Empty_Mangle()
    else:
        bleed_mask = pym.Mangle(config['fn_mask_bleed_pol'])

    ##### Check whether object are in the tolys:  since this is in a routine that will be call 5 times, it will be re-checked 5 times whereas once suffice ...

    tol = toly_mask.polyid(ra, dec)
    miscutils.fwdebug_print(f"There are {len(tol[tol < 0]):d} objects in the tiles which are lost by Mangle")

    ### Check whether objects are in the geometry mask.
    ## A[0] give the id of the polygon (-1 if no polygon), and A[1] give the associated mag_lim
    A = mangle_mask.polyid_and_weight(ra, dec)
    if miscutils.fwdebug_check(10, "MANGLE_DEBUG"):
        miscutils.fwdebug_print(f"A = {str(A)}")
        miscutils.fwdebug_print(f"A[1][A[0] < 10] = {A[1][A[0] < 10]}")

    star = star_mask.polyid(ra, dec)
    bleed = bleed_mask.polyid(ra, dec)

    if config['pltprefix'] is not None:
        plt.figure(1)
        plt.clf()
        plt.plot(ra[tol < 0], dec[tol < 0], 'k.')
        plt.savefig(f"{config['pltprefix']}_Object-not-in-tolygon.png")

        plt.figure(2)
        plt.clf()
        plt.plot(ra[A[0] < 0], dec[A[0] < 0], 'k.')
        plt.savefig(f"{config['pltprefix']}_Objects-not-in-geometry-masks.png")

        plt.figure(3)
        plt.clf()
        plt.plot(ra[star >= 0], dec[star >= 0], 'k.')
        plt.savefig(f"{config['pltprefix']}_Objects-in-star-mask.png")

        plt.figure(4)
        plt.clf()
        plt.plot(ra[bleed > 0], dec[bleed > 0], 'k.')
        plt.savefig(f"{config['pltprefix']}_Objects-in-bleed-mask.png")


    #### BIT values
    ISNOT_INGEOMETRY = 1
    IS_INSTAR = 2
    IS_INBLEED = 4
    ISNOT_INTOLY = 8

    fn_coadd_cat = config['fn_coadd_cat']
    if fn_coadd_cat is None:
        fn_coadd_cat = coadd_object_tab['FILENAME'][0]

    linefmt = '%d,%s,%s,%d,%d'
    g = open(fn_csv, 'w')
    for i in range(0, len(ra)):
        molygon_id_x = A[0][i]
        mangle_flag_x = ISNOT_INGEOMETRY * (A[0][i] < 0) + \
                        IS_INSTAR * (star[i] >= 0) + \
                        IS_INBLEED * (bleed[i] >= 0) + \
                        ISNOT_INTOLY * (tol[i] < 0)

        #print >>g, linefmt % (coadd_object_tab['ID'][i], molygon_id_x, mangle_flag_x)
        print(linefmt % (coadd_object_tab['NUMBER'][i], config['pfw_attempt_id'], config['band'],
                         molygon_id_x, mangle_flag_x), file=g)
    g.close()



######################################################################
def make_csv_ccdgon(config, fn_csv, Image_tab, log):

    coadd_filename = miscutils.parse_fullname(config['fn_coadd'], miscutils.CU_PARSE_FILENAME)

    jfn_weights = f"jcsvweights_{config['tileid']}_{config['band']}"
    cmd = f"poly2poly -ow {config['fn_ccdgon_pol']} {jfn_weights}"
    mu.runcmd(cmd, config['manglebindir'], log)
    weights = np.loadtxt(jfn_weights, skiprows=2)[:, 0]


    jfn_ids = f"jcsvids_{config['tileid']}_{config['band']}"
    cmd = f"poly2poly -oi {config['fn_ccdgon_pol']} {jfn_ids}"
    mu.runcmd(cmd, config['manglebindir'], log)
    if not os.path.exists(jfn_ids):
        raise IOError(f"Error could not find output file from poly2poly ({jfn_ids})")
    ids = np.loadtxt(jfn_ids, skiprows=2, dtype=np.int64)
    miscutils.fwdebug_print(f"Number of non-unique ccdgon polys = {len(ids)}")

    red_image_filename = Image_tab['FILENAME'][id2indices(Image_tab, ids)[0]]

    # eliminate ccdgons with duplicate ids
    _, idxs = np.unique(ids, return_index=True)
    miscutils.fwdebug_print(f"Number of unique ccdgon polys = {idxs.shape[0]}")

    id2idx = id2indices(Image_tab, ids)

    linefmt = '%d,%s,%s,%s,%f,%s,%s'
    g = open(fn_csv, 'w')
    for i in range(0, idxs.shape[0]):
        idx = idxs[i]
        ccdgon_id = ids[idx]
        print(i, idx, ccdgon_id)
        if id2idx[1][idx] == 0:
            amp = 'A'
        elif id2idx[1][idx] == 1:
            amp = 'B'
        else:
            raise KeyError(f"Invalid amp num for image #{idx} ({id2idx[1][idx]})")
        print(linefmt % (ccdgon_id, red_image_filename[idx], coadd_filename,
                         amp, 1./weights[idx], config['pfw_attempt_id'], config['band']), file=g)
    g.close()

    if config['cleanup'] is not None and config['cleanup'].upper() == 'Y':
        for tempf in [jfn_weights, jfn_ids]:
            if os.path.exists(tempf):
                os.remove(tempf)
            else:
                print(f"WARN: Could not find temporary file to delete ({tempf})")


######################################################################
def make_csv_molygon(config, fn_csv, molyids, log, Image_tab, Nmolys):

    jfn_mid = f"jweight-pol-mid_{config['tileid']}_{config['band']}"
    cmd = f"poly2poly -om {config['fn_molys_weight_pol']} {jfn_mid}"
    mu.runcmd(cmd, config['manglebindir'], log)
    midradec = np.loadtxt(jfn_mid, skiprows=2, usecols=[0, 1])

    ra_mid = midradec[:, 0]
    dec_mid = midradec[:, 1]
    coadd_filename = miscutils.parse_fullname(config['fn_coadd'], miscutils.CU_PARSE_FILENAME)
    coaddtile_id = config['tileid']

    band = config['band']

    make_syste_masks('EXPTIME', config, Nmolys, Image_tab, molyids)
    make_syste_masks('AIRMASS', config, Nmolys, Image_tab, molyids)
    make_syste_masks('SKYBRITE', config, Nmolys, Image_tab, molyids)
    make_syste_masks('SKYSIGMA', config, Nmolys, Image_tab, molyids)
    make_syste_masks('FWHM', config, Nmolys, Image_tab, molyids)
    make_syste_masks('T_EFF', config, Nmolys, Image_tab, molyids)
    make_syste_masks('FWHM_FLUXRAD', config, Nmolys, Image_tab, molyids)
    make_syste_masks('T_EFF_EXPTIME', config, Nmolys, Image_tab, molyids)
    make_syste_masks('SIGMA_MAG_ZERO', config, Nmolys, Image_tab, molyids, keys=['QSUM'])
    make_syste_masks('FGCM_GRY', config, Nmolys, Image_tab, molyids, exclude=0.0)
    make_syste_masks('SKYVAR', config, Nmolys, Image_tab, molyids, keys=['UNCERTAINTY', 'MIN', 'MAX'])

    num_images = np.loadtxt(config['fn_molys_count'])
    total_exptime = fitsio.read(config['data_file'], columns='SUM', ext='EXPTIME')
    area_str = np.loadtxt(config['fn_molys_area'])
    wavg_airmass = fitsio.read(config['data_file'], columns='WMEAN', ext='AIRMASS')
    wavg_fwhm = fitsio.read(config['data_file'], columns='WMEAN', ext='FWHM')
    mag_limit = np.loadtxt(config['fn_molys_maglims'])
    t_sigma_mag_zero = fitsio.read(config['data_file'], columns='QSUM', ext='SIGMA_MAG_ZERO')
    t_eff_exptime = fitsio.read(config['data_file'], columns='SUM', ext='T_EFF_EXPTIME')
    wavg_teff = fitsio.read(config['data_file'], columns='WMEAN', ext='T_EFF')
    wavg_fwhm_fluxr = fitsio.read(config['data_file'], columns='WMEAN', ext='FWHM_FLUXRAD')
    wavg_fgcm_gry = fitsio.read(config['data_file'], columns='WMEAN', ext='FGCM_GRY')
    min_fgcm_gry = fitsio.read(config['data_file'], columns='MIN', ext='FGCM_GRY')


    min_airmass = fitsio.read(config['data_file'], columns='MIN', ext='AIRMASS')
    max_airmass = fitsio.read(config['data_file'], columns='MAX', ext='AIRMASS')
    min_fwhm = fitsio.read(config['data_file'], columns='MIN', ext='FWHM')
    max_fwhm = fitsio.read(config['data_file'], columns='MAX', ext='FWHM')
    min_teff = fitsio.read(config['data_file'], columns='MIN', ext='T_EFF')
    max_teff = fitsio.read(config['data_file'], columns='MAX', ext='T_EFF')
    min_fwhm_fluxr = fitsio.read(config['data_file'], columns='MIN', ext='FWHM_FLUXRAD')
    max_fwhm_fluxr = fitsio.read(config['data_file'], columns='MAX', ext='FWHM_FLUXRAD')
    wavg_skyvar = fitsio.read(config['data_file'], columns='UNCERTAINTY', ext='SKYVAR')
    min_skyvar = fitsio.read(config['data_file'], columns='MIN', ext='SKYVAR')
    max_skyvar = fitsio.read(config['data_file'], columns='MAX', ext='SKYVAR')

    linefmt = '%d,%s,%s,%s,%d,%f,%15.7e,%f,%f,%f,%f,%f,%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f'

    if miscutils.fwdebug_check(6, "MANGLE_DEBUG"):
        i = 0
        miscutils.fwdebug_print(linefmt % (molyids[i], coadd_filename, coaddtile_id,
                                           band, num_images[i], total_exptime[i],
                                           area_str[i], wavg_airmass[i], wavg_fwhm[i],
                                           mag_limit[i], ra_mid[i], dec_mid[i],
                                           config['pfw_attempt_id'], t_sigma_mag_zero[i], t_eff_exptime[i], wavg_teff[i],
                                           wavg_fwhm_fluxr[i], wavg_fgcm_gry[i], min_fgcm_gry[i],
                                           min_airmass[i], max_airmass[i], min_fwhm[i], max_fwhm[i], min_teff[i],
                                           max_teff[i], min_fwhm_fluxr[i], max_fwhm_fluxr[i], wavg_skyvar[i],
                                           min_skyvar[i], max_skyvar[i]))

    g = open(fn_csv, 'w')
    for i in range(0, len(molyids)):
        #print i
        print(linefmt % (molyids[i], coadd_filename, coaddtile_id, band, num_images[i],
                         total_exptime[i], area_str[i], wavg_airmass[i], wavg_fwhm[i],
                         mag_limit[i], ra_mid[i], dec_mid[i], config['pfw_attempt_id'],
                         t_sigma_mag_zero[i], t_eff_exptime[i], wavg_teff[i],
                         wavg_fwhm_fluxr[i], wavg_fgcm_gry[i], min_fgcm_gry[i],
                         min_airmass[i], max_airmass[i], min_fwhm[i], max_fwhm[i], min_teff[i],
                         max_teff[i], min_fwhm_fluxr[i], max_fwhm_fluxr[i], wavg_skyvar[i],
                         min_skyvar[i], max_skyvar[i]), file=g)
    g.close()

    if config['cleanup'] is not None and config['cleanup'].upper() == 'Y':
        if os.path.exists(jfn_mid):
            os.remove(jfn_mid)
        else:
            print(f"WARN: Could not find temporary file to delete ({jfn_mid})")



######################################################################
def make_csv_molygon_ccdgon(config, fn_csv, molyids):

    linefmt = '%d,%d'

    g = open(fn_csv, 'w')
    redfh = open(config['fn_molys_red'])
    for i in range(0, len(molyids)):
        aa = redfh.readline().strip().split()
        n = len(aa) - 2
        for j in range(0, n):
            print(linefmt % (molyids[i], np.int(aa[2+j])), file=g)
    g.close()

######################################################################
def make_csv_files(config, Image_tab, dbi, Nmolys):
    """ Main function to create csv files for ingestion into db """

    if miscutils.fwdebug_check(10, "MANGLE_DEBUG"):
        miscutils.fwdebug_print(f"config = {config}")

    log = None
    if config['logdir'] is not None:
        log = config['logdir'] + 'make_csv_' + config['tilename'] + '_' + config['band'] +'.log'

    miscutils.fwdebug_print(f"Number of reduced images: {len(Image_tab['FILENAME'])}")

    ### get information about coadd objects from the DB
    coadd_object_tab = None
    if config['fn_coadd_cat'] is None:
        coadd_object_tab = mdb.get_coadd_object_info(dbi,
                                                     config['tilename'],
                                                     config['orig_pfw_id'],
                                                     config['ra_column'],
                                                     config['dec_column'],
                                                     config['schema'])
    else:
        coadd_object_tab = mu.read_coadd_object_cat(dbi,
                                                    config['fn_coadd_cat'],
                                                    config['ra_column'],
                                                    config['dec_column'])


    if coadd_object_tab is None or len(coadd_object_tab) == 0:
        raise Exception("Could not load data from coadd_object_tab.")
    fndb_pattern = config['dbprefix'] + '_%s.csv'

    #### maybe implement something to check whether there is a mangle mask in this band ?

    miscutils.fwdebug_print(f"Number of objects from coadd_object table: {len(coadd_object_tab[config['ra_column']])}")

    # get molygon ids
    jfn_ids = f"jweight-pol-ids_{config['tileid']}_{config['band']}"
    cmd = f"poly2poly -oi {config['fn_molys_weight_pol']} {jfn_ids}"
    mu.runcmd(cmd, config['manglebindir'], log)
    moly_ids = np.loadtxt(jfn_ids, skiprows=2, dtype=np.int64)
    miscutils.fwdebug_print(f"Number of ids from moly weight pol = {len(moly_ids)}")

    if config['cleanup'] is not None and config['cleanup'].upper() == 'Y':
        if os.path.exists(jfn_ids):
            os.remove(jfn_ids)
        else:
            print(f"WARN: Could not find temporary file to delete ({jfn_ids})")

    #################################
    #### make coadd_object_molygon .csv file
    #################################
    make_csv_coadd_object_molygon(config, fndb_pattern % 'coadd_object_molygon', coadd_object_tab)

    ##################################
    ##### Make CCDGON .csv file
    ###################################
    make_csv_ccdgon(config, fndb_pattern % 'ccdgon', Image_tab, log)

    ###########################
    ### Make MOLYGON .csv file
    ##########################
    make_csv_molygon(config, fndb_pattern % 'molygon', moly_ids, log, Image_tab, Nmolys)

    ###################################
    ### Make molygon_ccdgon .csv file
    ####################################
    make_csv_molygon_ccdgon(config, fndb_pattern % 'molygon_ccdgon', moly_ids)

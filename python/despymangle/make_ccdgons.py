
import numpy as np
import os
import re
import fitsio
import sys
from shutil import copyfile

import despydb.desdbi as desdbi
from despymisc import miscutils 
from despyastro import wcsutil

import despymangle.mangle_utils as mu
import despymangle.mangle_db as mdb

######################################################################
def split_datasec(dsec):
    """ Divide datasec string into individual numbers """

    # dsec = [1:1024,1:4096]
    secmatch = re.match(r'\[(\d+):(\d+)\s*,\s*(\d+):(\d+)\]', dsec)
    if not secmatch:
        raise ValueError("Invalid datasecA (%s)" % dsec)

    x0 = np.int(secmatch.group(1))
    x1 = np.int(secmatch.group(2))
    y0 = np.int(secmatch.group(3))
    y1 = np.int(secmatch.group(4))

    return (x0, x1, y0, y1)


######################################################################
def create_header_from_db_info(nwgint_tab, i):
    """ From nwgint info gathered from image table, create a header """

    raise Exception('Not implemented yet')

    # example code from Felipe
    #header = fitsio.FITSHDR()
    #for k, v in image_dict.items():
    #    new_record = {'name': k,'value':v[0]}
    #    header.add_record(new_record)
    #
    #return header



######################################################################
def ccdgons(config, Image_tab, nwgint_tab, head_tab, dbi=None):
    """ First pass to ??? """

    tilename = config['tilename']
    tileid = config['tileid']
    band = config['band']
    mzpglobal = config['mzpglobal']
    mtol = config['mtol']
    pfwidtemp = config['pfw_attempt_id']

    manglebindir = config['manglebindir']

    log = None
    if config['logdir'] is not None:
        log = config['logdir'] + 'make_ccdgons_' + tilename + '_' + band + '.log'

    ### The existance of those files have been checked previously
    tilefile = config['poltiles']

    #######  Get data for streaks.
    if dbi is None:
        dbi = desdbi.DesDbi(config['des_services'], config['db_section'], True)

    streak_tab = mdb.get_streak_info(dbi, config['schema'], Image_tab['FILENAME'])

    st_list = []  ## list of the individal ccdgons that are trimmed.
    N_images = len(Image_tab['FILENAME'])
    for i in range(N_images):   ## for each image in the image_tab_band
        print '########## working on image %s'%Image_tab['FILENAME'][i]

        if 'FULLNAME' in nwgint_tab: 
            h = np.where(np.logical_and((nwgint_tab['CCDNUM'] == Image_tab['CCDNUM'][i]), (nwgint_tab['EXPNUM'] == Image_tab['EXPNUM'][i])))
            if h == False:
                raise KeyError("Could not find matching nwgint for image %s" % Image_tab['FILENAME'][i])

            fits_file = nwgint_tab['FULLNAME'][h[0][0]]
            hdr = fitsio.read_header(fits_file)
        else:
            raise KeyError('Missing FULLNAME in nwgint_tab')
            #hdr = create_head_from_db_info(nwgint_tab, i)

        ### GET CCDNUM.
        CCDNUM = hdr['CCDNUM']


        b1 = int(config['border1'])
        b2 = int(config['border2'])

        # examples of datasec
        #DES2356+0043_r5p02_D00239652_i_c30_nwgint.fits
        #DATASECA= '[1025:2048,1:4096]' / Data section from amp A
        #DATASECB= '[1:1024,1:4096]'    / Data section from amp B
        #DES2356+0043_r5p02_D00239652_i_c31_nwgint.fits
        #DATASECA= '[1025:2048,1:4096]' / Data section from amp A
        #DATASECB= '[1:1024,1:4096]'    / Data section from amp B
        #DES2356+0043_r5p02_D00243480_z_c32_nwgint.fits
        #DATASECA= '[1:1024,1:4096]'    / Data section from amp A
        #DATASECB= '[1025:2048,1:4096]' / Data section from amp B
        if CCDNUM >= 32:
            dataseca='[1:1024,1:4096]'     # Data section from amp A
            datasecb='[1025:2048,1:4096]'  # Data section from amp B
        else:
            dataseca='[1025:2048,1:4096]' # Data section from amp A
            datasecb='[1:1024,1:4096]'     # Data section from amp B

        ### GET corners and middles RA, DEC.

        ##1) Do. Amp A
        if CCDNUM != 31:
            (x0,x1,y0,y1) = split_datasec(dataseca)

            if x0 == 1:
                xa = x0+b1
                ya = y0+b2
                xc = x1
                yc = y0+b2
                xe = x1
                ye = y1-b2
                xg = x0+b1
                yg = y1-b2

            if x0 == 1025:
                xa = x0
                ya = y0+b2
                xc = x1-b1
                yc = y0+b2
                xe = x1-b1
                ye = y1-b2
                xg = x0
                yg = y1-b2


            wcs = wcsutil.WCS(hdr)

            raA_a, decA_a = wcs.image2sky(xa, ya)
            raA_c, decA_c = wcs.image2sky(xc, yc)
            raA_e, decA_e = wcs.image2sky(xe, ye)
            raA_g, decA_g = wcs.image2sky(xg, yg)

            raA_b = mu.mean_ra(raA_a, raA_c)
            decA_b = np.mean([decA_a, decA_c])

            raA_d = mu.mean_ra(raA_c, raA_e)
            decA_d = np.mean([decA_c, decA_e])

            raA_f = mu.mean_ra(raA_e, raA_g)
            decA_f = np.mean([decA_e, decA_g])

            raA_h = mu.mean_ra(raA_g, raA_a)
            decA_h = np.mean([decA_g, decA_a])


        ##2) Do. Amp B
        (x0,x1,y0,y1) = split_datasec(datasecb)

        if x0 == 1:
            xa = x0+b1
            ya = y0+b2
            xc = x1
            yc = y0+b2
            xe = x1
            ye = y1-b2
            xg = x0+b1
            yg = y1-b2

        if x0 == 1025:
            xa = x0
            ya = y0+b2
            xc = x1-b1
            yc = y0+b2
            xe = x1-b1
            ye = y1-b2
            xg = x0
            yg = y1-b2

        wcs = wcsutil.WCS(hdr)

        raB_a, decB_a = wcs.image2sky(xa, ya)
        raB_c, decB_c = wcs.image2sky(xc, yc)
        raB_e, decB_e = wcs.image2sky(xe, ye)
        raB_g, decB_g = wcs.image2sky(xg, yg)

        raB_b = mu.mean_ra(raB_a, raB_c)
        decB_b = np.mean([decB_a, decB_c])

        raB_d = mu.mean_ra(raB_c, raB_e)
        decB_d = np.mean([decB_c, decB_e])

        raB_f = mu.mean_ra(raB_e, raB_g)
        decB_f = np.mean([decB_e, decB_g])

        raB_h = mu.mean_ra(raB_g, raB_a)
        decB_h = np.mean([decB_g, decB_a])

        ### check if there streaks here!
        if streak_tab == False:
            h = []
        else:
            h = np.where(streak_tab['REDFILENAME'] == Image_tab['FILENAME'][i])[0]

        if len(h) == 0: ## there is no streak for this image
            # Do ampA
            if CCDNUM != 31:
                jfn_edges0 = 'jedges_%s_%s_%d0' % (tileid, band, i)
                mu.create_weighted_poly_from_edges16(raA_a, decA_a, raA_b, decA_b,
                                                     raA_c, decA_c, raA_d, decA_d,
                                                     raA_e, decA_e, raA_f, decA_f,
                                                     raA_g, decA_g, raA_h, decA_h,
                                                     '%s0' % Image_tab['MANGLE_IMG_ID'][i],
                                                     hdr['SKYVARA']*100**((mzpglobal-Image_tab['MAG_ZERO'][i])/2.5),
                                                     jfn_edges0)
                if os.path.exists(jfn_edges0):
                    st_list.append(jfn_edges0)

            # Do ampB
            jfn_edges1 = 'jedges_%s_%s_%d1' % (tileid, band, i)
            mu.create_weighted_poly_from_edges16(raB_a, decB_a, raB_b, decB_b,
                                                 raB_c, decB_c, raB_d, decB_d,
                                                 raB_e, decB_e, raB_f, decB_f,
                                                 raB_g, decB_g, raB_h, decB_h,
                                                 '%s1' % Image_tab['MANGLE_IMG_ID'][i],
                                                 hdr['SKYVARB']*100**((mzpglobal-Image_tab['MAG_ZERO'][i])/2.5),
                                                 jfn_edges1)
            if os.path.exists(jfn_edges1):
                st_list.append(jfn_edges1)

        else: ### Streaks !!!
            print "Streaks !!!"
            if miscutils.fwdebug_check(3, "MAKECCDGONS_DEBUG"):
                miscutils.fwdebug_print("%s streaks for image %s (pos: %s)" % (len(h), Image_tab['FILENAME'][i], h))
            if miscutils.fwdebug_check(9, "MAKECCDGONS_DEBUG"):
                for nb_st in range(len(h)):
                    miscutils.fwdebug_print("streaks %s h[nb_st] = %s, tab=%s)" % (nb_st, h[nb_st], streak_tab[h[nb_st]]))

            # Do ampA
            jfn_edges0_t = 'jedges_%s_%s_%d_%d0' % (tileid, band, i, i)
            imageid = '%s0' % Image_tab['MANGLE_IMG_ID'][i]
            
            if CCDNUM != 31:
                mu.create_weighted_poly_from_edges16(raA_a, decA_a, raA_b, decA_b,
                                                     raA_c, decA_c, raA_d, decA_d,
                                                     raA_e, decA_e, raA_f, decA_f,
                                                     raA_g, decA_g, raA_h, decA_h,
                                                     imageid,
                                                     hdr['SKYVARA']*100**((mzpglobal-Image_tab['MAG_ZERO'][i])/2.5),
                                                     jfn_edges0_t)

            # Do ampB
            jfn_edges1_t = 'jedges_%s_%s_%d_%d1' % (tileid, band, i, i)
            mu.create_weighted_poly_from_edges16(raB_a, decB_a, raB_b, decB_b,
                                                 raB_c, decB_c, raB_d, decB_d,
                                                 raB_e, decB_e, raB_f, decB_f,
                                                 raB_g, decB_g, raB_h, decB_h,
                                                 '%s1' % Image_tab['MANGLE_IMG_ID'][i],
                                                 hdr['SKYVARB']*100**((mzpglobal-Image_tab['MAG_ZERO'][i])/2.5),
                                                 jfn_edges1_t)


            #### Defining the streaks
            stl = []
            for nb_st in range(len(h)):
                jfn_edges_st = 'jedges_st_%s_%s_%s_%d' % (tileid, band, i, nb_st)

                ## list of the streak polygons in case there are many of them
                if len(streak_tab.shape) > 0:
                    mu.create_weighted_poly_from_edges(streak_tab[h[nb_st]]['RA_1'], streak_tab[h[nb_st]]['DEC_1'],
                                                       streak_tab[h[nb_st]]['RA_4'], streak_tab[h[nb_st]]['DEC_4'],
                                                       streak_tab[h[nb_st]]['RA_3'], streak_tab[h[nb_st]]['DEC_3'],
                                                       streak_tab[h[nb_st]]['RA_2'], streak_tab[h[nb_st]]['DEC_2'],
                                                       999, -1, jfn_edges_st)

                    stl.append(jfn_edges_st)


            ##### Creating the polygon file for the ccdgond that has a streak ###
            jfn_pol = 'jpol_%s_%s_%s' % (tileid, band, i)
            if CCDNUM != 31:
                cmd = 'poly2poly %s %s %s %s' % (mtol, jfn_edges0_t, jfn_edges1_t, jfn_pol)
            else:
                cmd = 'poly2poly %s %s %s' % (mtol, jfn_edges1_t, jfn_pol)    #### This was the fix Aug. 4 2016. replaced jfn_edges1 by jfn_edges1_t
            mu.runcmd(cmd, manglebindir, log)

            jfn_ss = 'jss_%s_%s_%s' % (tileid, band, i)
            cmd = 'snap -S %s %s %s %s' % (config['snap'], mtol, jfn_pol, jfn_ss)
            mu.runcmd(cmd, manglebindir, log)

            ##### Creating the polygon file for the streak
            jfn_pol_st = 'jpol_st_%s_%s_%s' % (tileid, band, i)
            cmd = 'poly2poly %s %s %s' % (mtol, ' '.join(stl), jfn_pol_st)
            mu.runcmd(cmd, manglebindir, log)

            jfn_ss_st = 'jss_st_%s_%s_%s' % (tileid, band, i)
            cmd = 'snap -S %s %s %s %s' % (config['snap'], mtol, jfn_pol_st, jfn_ss_st)
            mu.runcmd(cmd, manglebindir, log)

            ### Balkanize the two file with -Bn to keep the lowest weight (-1) when there's a streak
            jfn_b = 'jb_%s_%s_%s' % (tileid, band, i)
            cmd = 'balkanize %s -Bn -vo %s %s %s' % (mtol, jfn_ss, jfn_ss_st, jfn_b)
            mu.runcmd(cmd, manglebindir, log)

            #### Poly2poly to remove the polygon with negative weight
            jfn_p = 'jp_%s_%s_%s' % (tileid, band, i)
            cmd = 'poly2poly %s -j0 %s %s' % (mtol, jfn_b, jfn_p)
            mu.runcmd(cmd, manglebindir, log)

            ##Snap as it seems it is always good to snap
            jfn_s = 'js_%s_%s_%s' % (tileid, band, i)
            cmd = 'snap -S %s %s %s' % (mtol, jfn_p, jfn_s)
            mu.runcmd(cmd, manglebindir, log)

            if os.path.exists(jfn_s):
                st_list.append(jfn_s)

            if config['cleanup'] is not None and config['cleanup'].upper() == 'Y':
                for tempf in stl:
                    if os.path.exists(tempf):
                        os.remove(tempf)
                    else:
                        print "WARN: Could not find temporary file to delete (%s)" % tempf
                if CCDNUM !=31:
                    for tempf in [jfn_pol, jfn_ss, jfn_pol_st, jfn_ss_st, 
                                  jfn_b, jfn_p, jfn_edges0_t, jfn_edges1_t]:
                        if os.path.exists(tempf):
                            os.remove(tempf)
                        else:
                            print "WARN: Could not find temporary file to delete (%s)" % tempf
                else:  ### Added this. Indeed it tried to remove some of the stuffs made for ccdnum-31 ampA whcih do not exist
                    for tempf in [jfn_pol, jfn_ss, jfn_pol_st, jfn_ss_st, 
                                  jfn_b, jfn_p, jfn_edges1_t]:
                        if os.path.exists(tempf):
                            os.remove(tempf)
                        else:
                            print "WARN: Could not find temporary file to delete (%s)" % tempf



    jfn_ss2 = 'jss2_%s_%s' % (tileid, band)
    cmd = 'snap -S %s %s %s %s' % (config['snap'], mtol, ' '.join(st_list), jfn_ss2)
    mu.runcmd(cmd, manglebindir, log)

    jfn_p = 'jp_%s_%s' % (tileid, band)
    cmd = 'pixelize -vo %s %s %s %s' % (config['pix'], mtol, jfn_ss2, jfn_p)
    mu.runcmd(cmd, manglebindir, log)

    jfn_s = 'js_%s_%s' % (tileid, band)
    cmd = 'snap %s %s %s %s' % (config['snap'], mtol, jfn_p, jfn_s)
    mu.runcmd(cmd, manglebindir, log)

    jfn_r = 'jr_%s_%s' % (tileid, band)
    cmd = 'rasterize -T %s %s %s %s' % (mtol, tilefile, jfn_s, jfn_r)
    mu.runcmd(cmd, manglebindir, log)
    
    copyfile(jfn_r, config['fn_ccdgon_pol'])

    if config['cleanup'] is not None and config['cleanup'].upper() == 'Y':
        for tempf in [jfn_ss2, jfn_p, jfn_s, jfn_r]:
            if os.path.exists(tempf):
                os.remove(tempf)
            else:
                print "WARN: Could not find temporary file to delete (%s)" % tempf
        for tempf in st_list:
            if os.path.exists(tempf):
                os.remove(tempf)
            else:
                print "WARN: Could not find temporary file to delete (%s)" % tempf

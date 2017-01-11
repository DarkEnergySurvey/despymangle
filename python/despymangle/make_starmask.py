import numpy as np
import os
import os.path
from shutil import copyfile

import despydb.desdbi as desdbi

import despymangle.mangle_utils as mu
import despymangle.mangle_db as mdb



######################################################################
def starmask(config, Image_tab, dbi):
    ### ABL : note, in this version the bleed mask and star mask are not balkanized. This saves a lot of time!


    # pulling some values from config as not to have to put config throughout code
    band = config['band']
    tilename = config['tilename']
    tilefile = config['poltiles']
    snap = config['snap']
    mtol = config['mtol']
    pix = config['pix']
    fn_mask_star = config['fn_mask_star_pol']
    fn_mask_bleed = config['fn_mask_bleed_pol']

    manglebindir = None
    if 'manglebindir' in config:
        manglebindir = config['manglebindir']

    log = None
    if config['logdir'] is not None:
        log = config['logdir'] + 'make_starmask_' + tilename + '_' + band +'.log'

    ### BLEED TRAIL

    if dbi is None:
        dbi = desdbi.DesDbi(config['des_services'], config['db_section'], True)

    #mdb.load_gtt_filename(dbi, Image_tab['FILENAME'])
    trailbox_tab = mdb.get_bleedtrail_info(dbi, config['schema'], Image_tab['FILENAME'])

    #N_trails = trailbox_tab.shape[0]
    N_images_with_trails = np.unique(trailbox_tab['FILENAME']).shape[0]

    bl_list = []
    for i in range(0, N_images_with_trails):
        N_t = np.max(trailbox_tab['RNUM'][np.where(trailbox_tab['FILENAME'] == np.unique(trailbox_tab['FILENAME'])[i])])
        for j in range(N_t):
            idx = np.where((trailbox_tab['RNUM'] == j+1)*  (trailbox_tab['FILENAME'] == np.unique(trailbox_tab['FILENAME'])[i]))[0][0]

            jfn_bleed_tr = 'jbleed_%s_%d_%d.pol' % (band, i, j)
            mu.create_weighted_poly_from_vertice(trailbox_tab['RA_1'][idx], trailbox_tab['DEC_1'][idx],
                                                 trailbox_tab['RA_2'][idx], trailbox_tab['DEC_2'][idx],
                                                 trailbox_tab['RA_3'][idx], trailbox_tab['DEC_3'][idx],
                                                 trailbox_tab['RA_4'][idx], trailbox_tab['DEC_4'][idx],
                                                 i, 1, jfn_bleed_tr)
            bl_list.append(jfn_bleed_tr)

    jfn_bleed1 = 'jbleed1%s.pol' % (band)
    cmd = 'snap -S %s %s %s %s' % (snap, mtol, ' '.join(bl_list), jfn_bleed1)
    mu.runcmd(cmd, manglebindir, log)

    jfn_bleed2 = 'jbleed2%s.pol' % (band)
    cmd = 'poly2poly -k1e-14,1e-2 %s %s' % (jfn_bleed1, jfn_bleed2)
    mu.runcmd(cmd, manglebindir, log)

    jfn_bleed3 = 'jbleed3%s.pol' % (band)
    cmd = 'pixelize %s %s %s %s' % (pix, mtol, jfn_bleed2, jfn_bleed3)
    mu.runcmd(cmd, manglebindir, log)

    jfn_bleed4 = 'jbleed4%s.pol' % (band)
    cmd = 'snap %s %s %s %s' % (snap, mtol, jfn_bleed3, jfn_bleed4)
    mu.runcmd(cmd, manglebindir, log)

    # remove tmp files
    if 'cleanup' in config and config['cleanup'].upper() == 'Y':
        for tempf in bl_list:
            if os.path.exists(tempf):
                os.remove(tempf)
            else:
                print "WARN: Could not find temporary file to delete (%s)" % tempf
        for tempf in [jfn_bleed1, jfn_bleed2, jfn_bleed3]:
            if os.path.exists(tempf):
                os.remove(tempf)
            else:
                print "WARN: Could not find temporary file to delete (%s)" % tempf

    jfn_rbleed = 'jrbleed%s.pol' % (band)
    cmd = 'rasterize -T %s %s %s %s' % (mtol, tilefile, jfn_bleed4, jfn_rbleed)
    mu.runcmd(cmd, manglebindir, log)

    # save under final filename
    copyfile(jfn_rbleed, fn_mask_bleed)

    # remove tmp files
    if 'cleanup' in config and config['cleanup'].upper() == 'Y':
        for tempf in [jfn_bleed4, jfn_rbleed]:
            if os.path.exists(tempf):
                os.remove(tempf)
            else:
                print "WARN: Could not find temporary file to delete (%s)" % tempf

    print  'bleedtrails done'

    #### STARS
    satstar_tab = mdb.get_satstar_info(dbi, config['schema'], Image_tab['FILENAME'], config['max_radius'])

    print 'je suis la'

    #N_stars = satstar_tab.shape[0]

    jfn_weightstar = 'jweightstar'
    os.system('echo 1 > %s' % jfn_weightstar)

    ufilenames = np.unique(satstar_tab['FILENAME'])
    N_images_with_stars = len(ufilenames)

    jfn_star = 'jstar_%s' % (band)
    fily = open(jfn_star, 'w')

    for i in range(0, N_images_with_stars):
        print i
        ## get number of stars in this image
        N_s = np.max(satstar_tab['RNUM'][np.where(satstar_tab['FILENAME'] == ufilenames[i])])

        for j in range(N_s):
            idx = np.where((satstar_tab['RNUM'] == j+1)*  (satstar_tab['FILENAME'] == ufilenames[i]))[0][0]
            print >>fily, satstar_tab['RA'][idx], satstar_tab['DEC'][idx], satstar_tab['RADIUS'][idx]/60./60. * config['starscale']
            ###        fily.close()

    fily.close()

    jfn_star1 = 'jstar1%s.pol' % (band)
    cmd = 'weight %s -ic1 -z%s %s %s' % (mtol, jfn_weightstar, jfn_star, jfn_star1)
    mu.runcmd(cmd, manglebindir, log)

    jfn_star2 = 'jstar2%s.pol' % (band)
    cmd = 'snap -S %s %s %s %s' % (snap, mtol, jfn_star1, jfn_star2)
    mu.runcmd(cmd, manglebindir, log)

    jfn_star3 = 'jstar3%s.pol' % (band)
    cmd = 'poly2poly -k1e-14,1e-2 %s %s' % (jfn_star2, jfn_star3)
    mu.runcmd(cmd, manglebindir, log)

    jfn_star4 = 'jstar4%s.pol' % (band)
    cmd = 'pixelize %s %s %s %s' % (pix, mtol, jfn_star3, jfn_star4)
    mu.runcmd(cmd, manglebindir, log)

    jfn_star5 = 'jstar5%s.pol' % (band)
    cmd = 'snap %s %s %s %s' % (snap, mtol, jfn_star4, jfn_star5)
    mu.runcmd(cmd, manglebindir, log)

    if 'cleanup' in config and config['cleanup'].upper() == 'Y':
        for tempf in [jfn_star, jfn_star1, jfn_star2, jfn_star3, jfn_star4]:
            if os.path.exists(tempf):
                os.remove(tempf)
            else:
                print "WARN: Could not find temporary file to delete (%s)" % tempf

    ###cmd = 'pixelize %s %s jstarss%s.pol jstarssp%s' % (pix, mtol, band, band)
    ###mu.runcmd(cmd, manglebindir, log)

    jfn_rstar = 'jrstar%s.pol' % (band)
    cmd = 'rasterize -T %s %s %s %s' % (mtol, tilefile, jfn_star5, jfn_rstar)
    mu.runcmd(cmd, manglebindir, log)

    # save data in final filename and delete intermediate files
    copyfile(jfn_rstar, fn_mask_star)
    if 'cleanup' in config and config['cleanup'].upper() == 'Y':
        for tempf in [jfn_star5, jfn_rstar, jfn_weightstar]:
            if os.path.exists(tempf):
                os.remove(tempf)
            else:
                print "WARN: Could not find temporary file to delete (%s)" % tempf


######################################################################
#def dectodeg(dec):
#    sign = np.sign(dec)
#    deg = np.floor(dec*sign)*sign
#    minu = np.floor(np.abs((dec-deg)) *60)
#    arcs =   np.abs((sign* dec-   sign*deg -  minu /60.))*3600
#    return deg, minu, arcs

######################################################################
#def distance(l1, l2):
#    return np.sqrt((l1[0]-l2[0])**2+  (l1[1]-l2[1])**2)

import os
from shutil import copyfile

import numpy as np

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

            jfn_bleed_tr = f"jbleed_{band}_{i:d}_{j:d}.pol"
            mu.create_weighted_poly_from_vertice(trailbox_tab['RA_1'][idx], trailbox_tab['DEC_1'][idx],
                                                 trailbox_tab['RA_2'][idx], trailbox_tab['DEC_2'][idx],
                                                 trailbox_tab['RA_3'][idx], trailbox_tab['DEC_3'][idx],
                                                 trailbox_tab['RA_4'][idx], trailbox_tab['DEC_4'][idx],
                                                 i, 1, jfn_bleed_tr)
            bl_list.append(jfn_bleed_tr)

    jfn_bleed1 = f"jbleed1{band}.pol"
    cmd = f"snap -S {snap} {mtol} {' '.join(bl_list)} {jfn_bleed1}"
    mu.runcmd(cmd, manglebindir, log)

    jfn_bleed2 = f"jbleed2{band}.pol"
    cmd = f"poly2poly -k1e-14,1e-2 {jfn_bleed1} {jfn_bleed2}"
    mu.runcmd(cmd, manglebindir, log)

    jfn_bleed3 = f"jbleed3{band}.pol"
    cmd = f"pixelize {pix} {mtol} {jfn_bleed2} {jfn_bleed3}"
    mu.runcmd(cmd, manglebindir, log)

    jfn_bleed4 = f"jbleed4{band}.pol"
    cmd = f"snap {snap} {mtol} {jfn_bleed3} {jfn_bleed4}"
    mu.runcmd(cmd, manglebindir, log)

    # remove tmp files
    if 'cleanup' in config and config['cleanup'].upper() == 'Y':
        for tempf in bl_list:
            if os.path.exists(tempf):
                os.remove(tempf)
            else:
                print(f"WARN: Could not find temporary file to delete ({tempf})")
        for tempf in [jfn_bleed1, jfn_bleed2, jfn_bleed3]:
            if os.path.exists(tempf):
                os.remove(tempf)
            else:
                print(f"WARN: Could not find temporary file to delete ({tempf})")

    jfn_rbleed = f"jrbleed{band}.pol"
    cmd = f"rasterize -T {mtol} {tilefile} {jfn_bleed4} {jfn_rbleed}"
    mu.runcmd(cmd, manglebindir, log)

    # save under final filename
    copyfile(jfn_rbleed, fn_mask_bleed)

    # remove tmp files
    if 'cleanup' in config and config['cleanup'].upper() == 'Y':
        for tempf in [jfn_bleed4, jfn_rbleed]:
            if os.path.exists(tempf):
                os.remove(tempf)
            else:
                print(f"WARN: Could not find temporary file to delete ({tempf})")

    print('bleedtrails done')

    #### STARS
    satstar_tab = mdb.get_satstar_info(dbi, config['schema'], Image_tab['FILENAME'], config['max_radius'])

    print('je suis la')

    #N_stars = satstar_tab.shape[0]

    jfn_weightstar = 'jweightstar'
    os.system(f"echo 1 > {jfn_weightstar}")

    ufilenames = np.unique(satstar_tab['FILENAME'])
    N_images_with_stars = len(ufilenames)

    jfn_star = f"jstar_{band}"
    fily = open(jfn_star, 'w')

    for i in range(0, N_images_with_stars):
        print(i)
        ## get number of stars in this image
        temp = np.where(satstar_tab['FILENAME'] == ufilenames[i])
        for idx in temp[0].tolist():
            print(satstar_tab['RA'][idx], satstar_tab['DEC'][idx], satstar_tab['RADIUS'][idx] / 60. / 60. * config['starscale'], file=fily)

    fily.close()

    if N_images_with_stars == 0:
        print("Found no satstars in database.")
        fl = open(fn_mask_star, 'w')
        fl.write('0 polygons\nreal 10\npixelization 10s\nsnapped\n')
        fl.close()
        return

    jfn_star1 = f"jstar1{band}.pol"
    cmd = f"weight {mtol} -ic1 -z{jfn_weightstar} {jfn_star} {jfn_star1}"
    mu.runcmd(cmd, manglebindir, log)

    jfn_star2 = f"jstar2{band}.pol"
    cmd = f"snap {snap} {mtol} {jfn_star1} {jfn_star2}"
    mu.runcmd(cmd, manglebindir, log)

    jfn_star3 = f"jstar3{band}.pol"
    cmd = f"poly2poly -k1e-14,1e-2 {jfn_star2} {jfn_star3}"
    mu.runcmd(cmd, manglebindir, log)

    jfn_star4 = f"jstar4{band}.pol"
    cmd = f"pixelize {pix} {mtol} {jfn_star3} {jfn_star4}"
    mu.runcmd(cmd, manglebindir, log)

    if 'cleanup' in config and config['cleanup'].upper() == 'Y':
        for tempf in [jfn_star, jfn_star1, jfn_star2, jfn_star3]:
            if os.path.exists(tempf):
                os.remove(tempf)
            else:
                print(f"WARN: Could not find temporary file to delete ({tempf})")

    ###cmd = 'pixelize %s %s jstarss%s.pol jstarssp%s' % (pix, mtol, band, band)
    ###mu.runcmd(cmd, manglebindir, log)

    jfn_rstar = f"jrstar{band}.pol"
    cmd = f"rasterize -T {mtol} {tilefile} {jfn_star4} {jfn_rstar}"
    mu.runcmd(cmd, manglebindir, log)

    # save data in final filename and delete intermediate files
    copyfile(jfn_rstar, fn_mask_star)
    if 'cleanup' in config and config['cleanup'].upper() == 'Y':
        for tempf in [jfn_star4, jfn_rstar, jfn_weightstar]:
            if os.path.exists(tempf):
                os.remove(tempf)
            else:
                print(f"WARN: Could not find temporary file to delete ({tempf})")


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

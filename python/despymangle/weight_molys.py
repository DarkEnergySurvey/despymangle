import numpy as np
import os
from shutil import copyfile

import despymangle.mangle_utils as mu

######################################################################
def id2skyvar(Image_tab, ids):
    ### return the bood skyvar whether the ccdgons is A or  B

    skyvars = []
    for one_id in ids:
        # drop last digit which means sec A vs B of image
        img_id = int(str(one_id)[:-1])

        is_AB=int(str(one_id)[-1])

        idx = np.where(Image_tab["MANGLE_IMG_ID"] == img_id)[0][0]
        if is_AB==0:
            skyvars.append(Image_tab["SKYVARA"][idx])
        if is_AB==1:
            skyvars.append(Image_tab["SKYVARB"][idx])
    return skyvars

######################################################################
def id2indices(Image_tab, ids):
    
    indices = []
    for one_id in ids:
        # drop last digit which means sec A vs B of image 
        img_id = int(str(one_id)[:-1])   

        # find the index in the Image_tab where the image with that mangle_img_id exists
        idx = np.where(Image_tab["MANGLE_IMG_ID"] == img_id)[0][0]
        indices.append(idx)

    return indices

######################################################################
def make_syste_masks(keyword, config, Nmolys, fn_reduced, Image_tab, fnprefix):
    # 6 as we want sum, min, max, median, mean and weightavg
    tab_exptime = np.zeros((6, Nmolys))
    f = open(fn_reduced)
    for i in range(Nmolys):   # for each molygon
        line = f.readline().strip()
        #ids = np.int32(np.array(line.split()[2:]))
        ids = line.split()[2:]    ### this all the ids of the ccdgons.
        indices = id2indices(Image_tab, ids)  ## indices is the indices in image_Tab which corresponds to the ccdgons
        skyvars= id2skyvar(Image_tab, ids)

        tab_exptime[0, i] = np.sum(Image_tab[keyword][indices])
        tab_exptime[1, i] = np.min(Image_tab[keyword][indices])
        tab_exptime[2, i] = np.max(Image_tab[keyword][indices])
        tab_exptime[3, i] = np.median(Image_tab[keyword][indices])
        tab_exptime[4, i] = np.mean(Image_tab[keyword][indices])

        weights = skyvars * 100 **   ((config['mzpglobal']-Image_tab['MAG_ZERO'][indices])/2.5)
        tab_exptime[5, i] = np.sum(Image_tab[keyword][indices] * weights)/np.sum(weights)

    f.close()
    fnprefix += '_' + keyword
    np.savetxt(fnprefix+'.SUM', tab_exptime[0, :])
    np.savetxt(fnprefix+'.MIN', tab_exptime[1, :])
    np.savetxt(fnprefix+'.MAX', tab_exptime[2, :])
    np.savetxt(fnprefix+'.MEDIAN', tab_exptime[3, :])
    np.savetxt(fnprefix+'.MEAN', tab_exptime[4, :])
    np.savetxt(fnprefix+'.WMEAN', tab_exptime[5, :])


######################################################################
def weightmolys(config, Image_tab):
    ###### stuff from weight files
    ###  In the new mangle pipeline, this can be done before starmask!!

    # pulling some values from config as not to have to put config throughout code
    band = config['band']
    tilename = config['tilename']
    tileid = config['tileid']
    #pfwidtemp = config['pfw_attempt_id']
    runn = config['runn']
    manglebindir = config['manglebindir']
    fn_maglimmask = config['fn_maglims']
    fn_weightmask = config['fn_molys_weight']

    log = None
    if 'logdir' in config and config['logdir'] is not None:
        log = config['logdir']+'weight_molys_'+tilename+'_'+band+'.log'

    print "calculating weights for molygons for %s tile %s (%s) %s band for run %s" % \
          (config['project'], tileid, tilename, band, runn)


    outdir = ''
    if config['outputdir'] is not None:
        outdir = config['outputdir'] + '/'
    fnprefix = outdir + config['molysprefix']
    fn_reduced = fnprefix + '.red'
    fn_count = fnprefix + '.count'
    fn_weight = fnprefix + '.weight'
    fn_maglims = fnprefix + '.maglims'
    fn_time = fnprefix + '.time'
    fn_area = fnprefix + '.area'

    #fn_bitmask = config['molysprefix'] + '.bitmask'
    #fn_starbitmask = config['outputdir']+'/%s_starbitmask_%s.pol'%(runn, band)
    ### copyfile(fn_starbitmask, 'jstarbitmask')

    #find midpoints of mask polygons

    jfn_mask = 'jmask_%s_%s' % (tileid, band)  ### jfn_mask is a copy of the ccdmolys_weight.pol
    copyfile(config['fn_mask'], jfn_mask)

    jfn_mid = 'jmid_%s_%s' % (tileid, band)    ### jfn_mid contains the midpoint of the molygons
    cmd = 'poly2poly -om %s %s' % (jfn_mask, jfn_mid)
    mu.runcmd(cmd, manglebindir, log)

    #use polyid to get id numbers of input reduced images
    jfn_ccdgon = 'jccdgon%s_%s' % (tileid, band)
    copyfile(config['fn_ccdgon'], jfn_ccdgon)

    jfn_f = 'jf_%s_%s' % (tileid, band)
    cmd = 'polyid %s %s %s' % (jfn_ccdgon, jfn_mid, jfn_f)
    mu.runcmd(cmd, manglebindir, log)

    jfn_reduced = 'jreduced_%s_%s' % (tileid, band)
    os.system('tail -n +2 %s > %s' % (jfn_f, jfn_reduced))

    #### get the nuumber of CCDgons per ccdmolygon
    jfn_count = 'jcount_%s_%s' % (tileid, band)
    os.system("awk '{print NF -2}' %s  > %s" % (jfn_reduced, jfn_count))

    copyfile(jfn_reduced, fn_reduced)
    copyfile(jfn_count, fn_count)

    if 'cleanup' in config and config['cleanup'].upper() == 'Y':
        for f in [jfn_reduced, jfn_count, jfn_f]:
            os.remove(f)

    #write areas
    jfn_a = 'ja_%s_%s' % (tileid, band)
    cmd = 'poly2poly -oa %s %s' % (jfn_mask, jfn_a)
    mu.runcmd(cmd, manglebindir, log)

    jfn_area = 'jarea_%s_%s' % (tileid, band)
    os.system("awk 'NR>2{print $1}' %s > %s " % (jfn_a, jfn_area))
    copyfile(jfn_area, fn_area)
    if 'cleanup' in config and config['cleanup'].upper() == 'Y':
        for junkf in [jfn_area, jfn_a]:
            os.remove(junkf)

    jfn_w = 'jw_%s_%s' % (tileid, band)

    ### load jcount
    jcount = np.loadtxt(fn_count)
    Nmolys = len(jcount)

### Do coadd weighting

    ### get for each center of molygons the lists of SE weights
    jfn_f = 'jSE_w_%s_%s' % (tileid, band)
    cmd = 'polyid -W %s %s %s' % (jfn_ccdgon, jfn_mid, jfn_f)
    mu.runcmd(cmd, manglebindir, log)

    jfn_fs = 'jfs_%s_%s' % (tileid, band)
    os.system('tail -n +2 %s > %s' % (jfn_f, jfn_fs))

    aaa = open(jfn_fs).readlines()
    Nmolys = len(aaa)
    print  'There are %d molys in this tile'%Nmolys
    weight_tot = np.zeros(Nmolys)

    ####abl  mzp = 30.0    

    for i in range(Nmolys):
        linee = aaa[i].strip().split()
        bbb = [1.0/np.float128(linee[j]) for j in range(2, len(linee))]
        ###abl  zzz[i] = 100**((config['mzpglobal']-np.float128(mzp))/2.5)*np.sum(bbb)
        weight_tot[i] = np.sum(bbb)

    jfn_weight = 'jweight_%s_%s' % (tileid, band)
    np.savetxt(jfn_weight, weight_tot, fmt='%.18e')
    cmd = 'weight -z%s %s %s' % (jfn_weight, jfn_mask, fn_weightmask)
    mu.runcmd(cmd, manglebindir, log)
    copyfile(jfn_weight, fn_weight)

    print  "wrote weights to %s" % fn_weight

    #calculate magnitude limit from coadded weight

    maglims = np.float128(config['mzpglobal']) - \
          2.5*np.log10(10*np.sqrt(np.pi)*(config['aper']/2/config['asperpix'])) - \
          2.5*np.log10(1.0/np.sqrt(weight_tot))

    jfn_maglims = 'jmaglims_%s_%s' % (tileid, band)
    np.savetxt(jfn_maglims, maglims, fmt='%.18e')
    copyfile(jfn_maglims, fn_maglims)

    print "wrote magnitude limits to %s" % fn_maglims

    jfn_jmaglimmask = 'jmaglimmask_%s_%s' % (tileid, band)
    cmd = 'weight -z%s %s %s' % (jfn_maglims, jfn_mask, jfn_jmaglimmask)
    mu.runcmd(cmd, manglebindir, log)
    copyfile(jfn_jmaglimmask, fn_maglimmask)

    make_syste_masks('EXPTIME', config, Nmolys, fn_reduced, Image_tab, fnprefix)
    make_syste_masks('AIRMASS', config, Nmolys, fn_reduced, Image_tab, fnprefix)
    make_syste_masks('SKYBRITE', config, Nmolys, fn_reduced, Image_tab, fnprefix)
    make_syste_masks('SKYSIGMA', config, Nmolys, fn_reduced, Image_tab, fnprefix)
    make_syste_masks('FWHM', config, Nmolys, fn_reduced, Image_tab, fnprefix)
    #make_syste_masks('SKYVARA', config, Nmolys, fn_reduced, Image_tab, fnprefix)
    #make_syste_masks('SKYVARB', config, Nmolys, fn_reduced, Image_tab, fnprefix)


    #write total observation times

    cmd = 'poly2poly -ow %s %s'%(jfn_mask, jfn_w)
    mu.runcmd(cmd, manglebindir, log)

    jfn_time = 'jtime_%s_%s' % (tileid, band) 
    os.system("awk 'NR>2{print $1}' %s > %s" % (jfn_w, jfn_time))
    copyfile(jfn_time, fn_time)

    if 'cleanup' in config and config['cleanup'].upper() == 'Y':
        for junkf in [jfn_time, jfn_w]:
            os.remove(junkf)

### Do coadd weigthing
###    cmd = 'polyid -W %s %s %s' % (jfn_ccdgon, jfn_mid, jfn_f)
###    mu.runcmd(cmd, manglebindir, log)
###    os.system('tail -n +2 %s >  %s' % (jfn_f, jfn_fs))

    # clean up temp files
    if 'cleanup' in config and config['cleanup'].upper() == 'Y':
        for tempf in [jfn_ccdgon, jfn_jmaglimmask, jfn_weight,
                      jfn_mask, jfn_maglims,
                      jfn_f, jfn_fs, jfn_mid]:
            os.remove(tempf)

import numpy as np
import os
from shutil import copyfile

import despymangle.mangle_utils as mu

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
def make_syste_masks(keyword, Nmolys, fn_reduced, Image_tab, fnprefix):
    # 6 as we want sum, min, max, median, mean and weightavg
    tab_exptime = np.zeros((6, Nmolys))
    f = open(fn_reduced)
    for i in range(Nmolys):
        line = f.readline().strip()
        #ids = np.int32(np.array(line.split()[2:]))
        ids = line.split()[2:]
        indices = id2indices(Image_tab, ids)
        tab_exptime[0, i] = np.sum(Image_tab[keyword][indices])
        tab_exptime[1, i] = np.min(Image_tab[keyword][indices])
        tab_exptime[2, i] = np.max(Image_tab[keyword][indices])
        tab_exptime[3, i] = np.median(Image_tab[keyword][indices])
        tab_exptime[4, i] = np.mean(Image_tab[keyword][indices])
        ####       tab_exptime[5, i] = np.sum(Image_tab[keyword][indices])

    f.close()
    fnprefix += '_' + keyword
    np.savetxt(fnprefix+'.SUM', tab_exptime[0, :])
    np.savetxt(fnprefix+'.MIN', tab_exptime[1, :])
    np.savetxt(fnprefix+'.MAX', tab_exptime[2, :])
    np.savetxt(fnprefix+'.MEDIAN', tab_exptime[3, :])
    np.savetxt(fnprefix+'.MEAN', tab_exptime[4, :])
    ##np.savetxt(fnprefix+'.WMEAN', tab_exptime[5, :])


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

    log = None
    if 'logdir' in config and config['logdir'] is not None:
        #log = logdir+'weight_molys_'+runn+'_'+band+'.log'
        # why was this one using runn instead of tilename like others?
        log = config['logdir']+'weight_molys_'+tilename+'_'+band+'.log'

    print "calculating weights for molygons for %s tile %s (%s) %s band for run %s" % \
          (config['project'], tileid, tilename, band, runn)

    fn_maglimmask = config['outputdir']+'/%s_molys_maglims_%s.pol'%(runn, band)
    fn_weightmask = config['outputdir']+'/%s_molys_weight_%s.pol'%(runn, band)

    fn_reduced = config['outputdir']+'/%s_molys_%s.red'%(runn, band)
    fn_count = config['outputdir']+'/%s_molys_%s.count'%(runn, band)
    fn_weight = config['outputdir']+'/%s_molys_%s.weight'%(runn, band)
    fn_maglims = config['outputdir']+'/%s_molys_%s.maglims'%(runn, band)
    fn_time = config['outputdir']+'/%s_molys_%s.time'%(runn, band)
    fn_area = config['outputdir']+'/%s_molys_%s.area'%(runn, band)

    #fn_bitmask = config['outputdir']+'/%s_molys_%s.bitmask'%(runn, band)
    #fn_starbitmask = config['outputdir']+'/%s_starbitmask_%s.pol'%(runn, band)
    ### copyfile(fn_starbitmask, 'jstarbitmask')

    #find midpoints of mask polygons

    jfn_mask = 'jmask_%s_%s' % (tileid, band)
    copyfile(config['fn_mask'], jfn_mask)

    jfn_mid = 'jmid_%s_%s' % (tileid, band)
    cmd = 'poly2poly -om %s %s' % (jfn_mask, jfn_mid)
    mu.runcmd(cmd, manglebindir, log)

    #use polyid to get id numbers of input reduced images
    jfn_unbalk = 'junbalk_%s_%s' % (tileid, band)
    copyfile(config['fn_unbalk'], jfn_unbalk)

    jfn_f = 'jf_%s_%s' % (tileid, band)
    cmd = 'polyid %s %s %s' % (jfn_unbalk, jfn_mid, jfn_f)
    mu.runcmd(cmd, manglebindir, log)

    jfn_reduced = 'jreduced_%s_%s' % (tileid, band)
    os.system('tail -n +2 %s > %s' % (jfn_f, jfn_reduced))

    jfn_count = 'jcount_%s_%s' % (tileid, band)
    os.system("awk '{print NF}' %s  > %s" % (jfn_reduced, jfn_count))

    os.rename(jfn_reduced, fn_reduced)
    os.rename(jfn_count, fn_count)

    if 'cleanup' in config and config['cleanup'].upper() == 'Y':
        os.remove(jfn_f)

    #write areas
    jfn_a = 'ja_%s_%s' % (tileid, band)
    cmd = 'poly2poly -oa %s %s' % (jfn_mask, jfn_a)
    mu.runcmd(cmd, manglebindir, log)

    jfn_area = 'jarea_%s_%s' % (tileid, band)
    os.system("awk 'NR>2{print $1}' %s > %s " % (jfn_a, jfn_area))
    os.rename(jfn_area, fn_area)
    if 'cleanup' in config and config['cleanup'].upper() == 'Y':
        os.remove(jfn_a)

    jfn_w = 'jw_%s_%s' % (tileid, band)

    ### load jcount
    jcount = np.loadtxt(fn_count)
    Nmolys = len(jcount)

### Do coadd weigthing

    ### get for each center of molygons the lists of SE weights
    cmd = 'polyid -W %s %s %s' % (jfn_unbalk, jfn_mid, jfn_f)
    mu.runcmd(cmd, manglebindir, log)

    jfn_fs = 'jfs_%s_%s' % (tileid, band)
    os.system('tail -n +2 %s > %s' % (jfn_f, jfn_fs))

    aaa = open(jfn_fs).readlines()
    Nmolys = len(aaa)
    print  'There are %d molys in this tile'%Nmolys
    zzz = np.zeros(Nmolys)

    mzp = 30.0    #MMG Should this be a value passed in - mag_zero or zp for each image?

    for i in range(Nmolys):
        linee = aaa[i].strip(). split()
        bbb = [1.0/np.float128(linee[j]) for j in range(2, len(linee))]
        zzz[i] = 100**((config['mzpglobal']-np.float128(mzp))/2.5)*np.sum(bbb)

    jfn_weight = 'jweight_%s_%s' % (tileid, band)
    np.savetxt(jfn_weight, zzz, fmt='%.18e')
    cmd = 'weight -z%s %s %s' % (jfn_weight, jfn_mask, fn_weightmask)
    mu.runcmd(cmd, manglebindir, log)
    copyfile(jfn_weight, fn_weight)

    print  "wrote weights to %s" % fn_weight

    #calculate magnitude limit from coadded weight

    yyy = np.float128(mzp) - \
          2.5*np.log10(10*np.sqrt(np.pi)*(config['aper']/2/config['asperpix'])) - \
          2.5*np.log10(1.0/np.sqrt(zzz))

    jfn_maglims = 'jmaglims_%s_%s' % (tileid, band)
    np.savetxt(jfn_maglims, yyy, fmt='%.18e')
    copyfile(jfn_maglims, fn_maglims)

    print "wrote magnitude limits to %s" % fn_maglims

    jfn_jmaglimmask = 'jmaglimmask_%s_%s' % (tileid, band)
    cmd = 'weight -z%s %s %s' % (jfn_maglims, jfn_mask, jfn_jmaglimmask)
    mu.runcmd(cmd, manglebindir, log)
    copyfile(jfn_jmaglimmask, fn_maglimmask)



    fnprefix = config['outputdir']+'/%s_molys_%s' % (runn, band)
    make_syste_masks('EXPTIME', Nmolys, fn_reduced, Image_tab, fnprefix)
    make_syste_masks('AIRMASS', Nmolys, fn_reduced, Image_tab, fnprefix)
    make_syste_masks('SKYBRITE', Nmolys, fn_reduced, Image_tab, fnprefix)
    make_syste_masks('SKYSIGMA', Nmolys, fn_reduced, Image_tab, fnprefix)
    make_syste_masks('FWHM', Nmolys, fn_reduced, Image_tab, fnprefix)
    make_syste_masks('SKYVARA', Nmolys, fn_reduced, Image_tab, fnprefix)
    make_syste_masks('SKYVARB', Nmolys, fn_reduced, Image_tab, fnprefix)


    #write total observation times

    cmd = 'poly2poly -ow %s %s'%(jfn_mask, jfn_w)
    mu.runcmd(cmd, manglebindir, log)

    os.system("awk 'NR>2{print $1}' %s > %s" % (jfn_w, 'jtime'))
    os.rename('jtime', fn_time)
    if 'cleanup' in config and config['cleanup'].upper() == 'Y':
        os.remove(jfn_w)

### Do coadd weigthing
    cmd = 'polyid -W %s %s %s' % (jfn_unbalk, jfn_mid, jfn_f)
    mu.runcmd(cmd, manglebindir, log)
    os.system('tail -n +2 %s >  %s' % (jfn_f, jfn_fs))

    aaa = open(jfn_fs).readlines()
    Nmolys = len(aaa)
    print  'There are %d molys in this tile'%Nmolys
    zzz = np.zeros(Nmolys)

    mzp = 30.0   #MMG already defined above?

    for i in range(Nmolys):
        linee = aaa[i].strip(). split()
        bbb = [1.0/np.float128(linee[j]) for j in range(2, len(linee))]
        zzz[i] = 100**((config['mzpglobal']-np.float128(mzp))/2.5)*np.sum(bbb)

    np.savetxt(jfn_weight, zzz, fmt='%.18e')
    cmd = 'weight -z%s %s %s' % (jfn_weight, jfn_mask, fn_weightmask)
    mu.runcmd(cmd, manglebindir, log)
    copyfile(jfn_weight, fn_weight)

    print  "wrote weights to %s" % fn_weight

    #calculate magnitude limit from coadded weight
    yyy = np.float128(mzp) - \
          2.5*np.log10(10*np.sqrt(np.pi)*(config['aper']/2/config['asperpix'])) - \
          2.5*np.log10(1.0/np.sqrt(zzz))

    np.savetxt(jfn_maglims, yyy, fmt='%.18e')
    copyfile(jfn_maglims, fn_maglims)

    print "wrote magnitude limits to %s" % fn_maglims

    cmd = 'weight -z%s %s %s' % (jfn_maglims, jfn_mask, jfn_jmaglimmask)
    mu.runcmd(cmd, manglebindir, log)
    copyfile(jfn_jmaglimmask, fn_maglimmask)

    # clean up temp files
    if 'cleanup' in config and config['cleanup'].upper() == 'Y':
        for tempf in [jfn_unbalk, jfn_jmaglimmask, jfn_weight,
                      jfn_mask, jfn_maglims,
                      jfn_f, jfn_fs, jfn_mid]:
            os.remove(tempf)

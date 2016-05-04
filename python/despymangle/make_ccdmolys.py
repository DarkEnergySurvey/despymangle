import os
from shutil import copyfile

import despymangle.mangle_utils as mu


def ccdmolys(config):

    # pulling some values from config as not to have to put config throughout code
    band = config['band']
    tilename = config['tilename']
    tileid = config['tileid']
    pfwidtemp = config['pfw_attempt_id']
    mtol = config['mtol']

    log = None
    if config['logdir'] is not None:
        log = config['logdir']+'make_ccdmolys_'+tilename+'_'+band+'.log'

    manglebindir = config['manglebindir']

    jfn_unbalk = 'junbalk_%s_%s' % (tileid, band)
    copyfile(config['fn_unbalk'], jfn_unbalk)

    jfn_b = 'jb_%s_%s' % (tileid, band)
    cmd = 'balkanize -Ba %s %s %s' % (mtol, jfn_unbalk, jfn_b)
    mu.runcmd(cmd, manglebindir, log)

    jfn_u = 'ju_%s_%s' % (tileid, band)
    cmd = 'unify %s %s' % (jfn_b, jfn_u)
    mu.runcmd(cmd, manglebindir, log)

    # get rid of polygons with areas less than 5e-14 str or greater than 6 str
    jfn_mask = 'jmask_%s_%s' % (tileid, band)
    numscheme = '%s%s000000' % (pfwidtemp, config['bandnum'])
    cmd = 'poly2poly -vn%s -k5e-14,6 %s %s' % (numscheme, jfn_u, jfn_mask)
    mu.runcmd(cmd, manglebindir, log)

    os.rename(jfn_mask, config['fn_mask'])

    if config['cleanup'] is not None and config['cleanup'].upper() == 'Y':
        for tempf in [jfn_b, jfn_u, jfn_unbalk]:
            os.remove(tempf)

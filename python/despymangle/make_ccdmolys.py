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
        log = config['logdir'] + 'make_ccdmolys_' + tilename + '_' + band + '.log'

    manglebindir = config['manglebindir']

    jfn_ccdgon = f"jccdgon{tileid}_{band}"
    copyfile(config['fn_ccdgon_pol'], jfn_ccdgon)

    jfn_b = f"jb_{tileid}_{band}"
    cmd = f"balkanize -Ba {mtol} {jfn_ccdgon} {jfn_b}"
    mu.runcmd(cmd, manglebindir, log)

    jfn_u = f"ju_{tileid}_{band}"
    cmd = f"unify {jfn_b} {jfn_u}"
    mu.runcmd(cmd, manglebindir, log)

    # get rid of polygons with areas less than 5e-14 str or greater than 6 str
    jfn_mask = f"jmask_{tileid}_{band}"
    numscheme = f"{pfwidtemp}{config['bandnum']}00000000"
    cmd = f"poly2poly -vn{numscheme} -k5e-14,6 {jfn_u} {jfn_mask}"
    mu.runcmd(cmd, manglebindir, log)

    copyfile(jfn_mask, config['fn_mask_pol'])

    if config['cleanup'] is not None and config['cleanup'].upper() == 'Y':
        for tempf in [jfn_b, jfn_u, jfn_ccdgon, jfn_mask]:
            if os.path.exists(tempf):
                os.remove(tempf)
            else:
                print(f"WARN: Could not find temporary file to delete ({tempf})")

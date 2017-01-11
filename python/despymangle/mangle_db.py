#! /usr/bin/env python

""" DB functions used by mangle """

import despyastro
from despymisc import miscutils
import despydb


###################################################################
def load_gtt_filename(dbi, redfiles):
    """ load ids into opm_filename_gtt table so can be used in sql joins """

    curs = dbi.cursor()

    # make sure opm_filename_gtt table is empty
    sql = "delete from OPM_FILENAME_GTT"
    curs = dbi.cursor()
    curs.execute(sql)

    # ensure only have filename (as DB column defined) instead of fullnames
    fnlist2 = []
    for fname in redfiles:
        fnlist2.append([miscutils.parse_fullname(fname, miscutils.CU_PARSE_FILENAME)])

    # load img ids into opm_filename_gtt table
    dbi.insert_many('OPM_FILENAME_GTT', ['FILENAME'], fnlist2)

###################################################################
def load_gtt_red_from_nwgint(dbi, nwgintfiles=None):
    """ Read nwgint information from database """

    # first load nwgint filenames into gtt 
    # in order to find red filenames
    load_gtt_filename(dbi, nwgintfiles)

    sql = """select dr.filename from desfile dr, desfile dn, 
             opm_was_derived_from od, opm_filename_gtt g 
             where dn.filename=g.filename and od.child_desfile_id=dn.id and 
             dr.id=parent_desfile_id and dr.filetype='red_immask'
          """ 

    curs = dbi.cursor()
    curs.execute(sql)
    redfiles = []
    for row in curs:
        redfiles.append(row[0])

    # now load red filenames into gtt for future queries
    load_gtt_filename(dbi, redfiles)
   

###################################################################
def get_redimg_info(dbi, schema='', redfiles=None):
    """ Read image information from database for reduced images """

    if redfiles is not None:
        load_gtt_filename(dbi, redfiles)
    sql = """select i.*, pq.fwhm_mean*0.263 as FWHM from image i, opm_filename_gtt g,miscfile m, 
             psf_qa pq where i.filename=g.filename and m.pfw_attempt_id=i.pfw_attempt_id and
             m.ccdnum=i.ccdnum and m.filename=pq.filename and m.filetype='xml_psfex' """

    if miscutils.fwdebug_check(3, "MANGLEDB_DEBUG"):
        miscutils.fwdebug_print("sql = %s" % sql)

    red_tab = despyastro.query2dict_of_columns(sql, dbi, True)
    return red_tab


###################################################################
def get_streak_info(dbi, schema='', redfiles=None):
    """ Read streak information from database """

    if redfiles is not None:
        load_gtt_filename(dbi, redfiles)
    
    sql = """select i.filename as redfilename, s.*
             from streak s, catalog c, image i, opm_filename_gtt g where
             s.filename=c.filename and i.filename=g.filename and
             i.expnum=c.expnum and i.ccdnum=c.ccdnum and i.pfw_attempt_id=c.pfw_attempt_id
          """

    if miscutils.fwdebug_check(3, "MANGLEDB_DEBUG"):
        miscutils.fwdebug_print("sql = %s" % sql)

    streak_tab = despyastro.query2rec(sql, dbhandle=dbi)
    return streak_tab


###################################################################
def get_satstar_info(dbi, schema='', redfiles=None, max_radius=400.0):
    """ Read satstar information from database """

    if redfiles is not None:
        load_gtt_filename(dbi, redfiles)
    
    sql = """select i.filename as redfilename, s.*
             from satstar s, catalog c, image i, opm_filename_gtt g where
             s.filename=c.filename and s.radius<%f and i.filename=g.filename and
             i.expnum=c.expnum and i.ccdnum=c.ccdnum and i.pfw_attempt_id=c.pfw_attempt_id
          """ % (max_radius)

    if miscutils.fwdebug_check(3, "MANGLEDB_DEBUG"):
        miscutils.fwdebug_print("sql = %s" % sql)

    satstar_tab = despyastro.query2rec(sql, dbhandle=dbi)
    return satstar_tab


###################################################################
def get_bleedtrail_info(dbi, schema='', redfiles=None):
    """ Read bleedtrail information from database """

    if redfiles is not None:
        load_gtt_filename(dbi, redfiles)
    
    sql = """select i.filename as redfilename, b.*
             from bleedtrail b, catalog c, image i, opm_filename_gtt g where
             b.filename=c.filename and i.filename=g.filename and
             i.expnum=c.expnum and i.ccdnum=c.ccdnum and i.pfw_attempt_id=c.pfw_attempt_id
          """

    if miscutils.fwdebug_check(3, "MANGLEDB_DEBUG"):
        miscutils.fwdebug_print("sql = %s" % sql)

    bleedtrail_tab = despyastro.query2rec(sql, dbhandle=dbi)
    return bleedtrail_tab

###################################################################
def get_headfile_info(dbi, me_attid, redfiles=None):
    """ Read coadd head file information from database """

    if redfiles is not None:
        load_gtt_filename(dbi, redfiles)

    sql = """select i.filename as redfilename, m.filename
             from miscfile m, image i, opm_filename_gtt g where
             i.filename=g.filename and
             i.expnum=m.expnum and i.ccdnum=m.ccdnum and m.pfw_attempt_id=%s
             and m.filetype='coadd_head_scamp'
          """ % (me_attid)

    if miscutils.fwdebug_check(3, "MANGLEDB_DEBUG"):
        miscutils.fwdebug_print("sql = %s" % sql)

    headfile_tab = despyastro.query2dict_of_columns(sql, dbi, True)
    return headfile_tab

##################################################################
def get_nwgint_info_from_reds(dbi, me_attid, redfiles=None):
    """ Read nwgint information from database """

    if redfiles is not None:
        load_gtt_filename(dbi, redfiles)
    
    sql = """select i1.filename as redfilename, i2.*
             from image i1, image i2, opm_filename_gtt g where
             i1.filename=g.filename and
             i1.expnum=i2.expnum and i1.ccdnum=i2.ccdnum and i2.pfw_attempt_id=%s
             and i2.filetype='coadd_nwgint'
          """ % (me_attid)

    if miscutils.fwdebug_check(3, "MANGLEDB_DEBUG"):
        miscutils.fwdebug_print("sql = %s" % sql)

    nwgint_tab = despyastro.query2dict_of_columns(sql, dbi, True)
    return nwgint_tab


###################################################################
def get_nwgint_info(dbi, schema='', nwgintfiles=None):
    """ Read nwgint information from database """

    if nwgintfiles is not None:
        load_gtt_filename(dbi, nwgintfiles)
    sql = """select i.* from image i, opm_filename_gtt g where 
             i.filename=g.filename and i.filetype='coadd_nwgint'
          """ 

    if miscutils.fwdebug_check(3, "MANGLEDB_DEBUG"):
        miscutils.fwdebug_print("sql = %s" % sql)

    nwgint_tab = despyastro.query2dict_of_columns(sql, dbi, True)
    return nwgint_tab


###################################################################
def get_tileid(tilename, db_section=None, des_services=None):
    """ Get tile id from database for given tilename """

    dbh = despydb.desdbi.DesDbi(des_services, db_section, True)

    # Prepare the query
    query = """select id from coaddtile_geom where tilename='%s'""" % tilename

    # Get the cursor
    cur = dbh.cursor()
    cur.execute(query)

    item = cur.fetchone()
    return item[0]


###################################################################
def get_coadd_object_info(dbi, tilename, pfw_attempt_id,
                          ra_column, dec_column, dbschema=''):
    """ Get id and location information for coadd objects """

    
    sql = """select co.id, co.filename, co.object_number as 'number', 
                    co.{racol} as ra, co.{deccol} as dec 
             from {schema}coadd_object co, {schema}catalog c 
             where co.filename=c.filename and 
                   c.filetype='coadd_det_cat' and 
                   c.tilename='{tname}' and 
                   c.pfw_attempt_id={attid}
          """.format(racol=ra_column, 
                     deccol=dec_column, 
                     schema=dbschema, 
                     tname=tilename, 
                     attid=pfw_attempt_id)

    if miscutils.fwdebug_check(3, "MANGLEDB_DEBUG"):
        miscutils.fwdebug_print("sql = %s" % sql)

    coadd_object_tab = despyastro.query2rec(sql, dbhandle=dbi)
    return coadd_object_tab


###################################################################
if __name__ == '__main__':
    pass

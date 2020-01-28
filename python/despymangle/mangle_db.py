""" DB functions used by mangle """

import numpy
import despyastro
from despymisc import miscutils
import despydb

###################################################################
def query2rec(query, dbhandle, verb=False):
    """
        Queries DB and returns results as a numpy recarray.
    """

    # Get the cursor from the DB handle
    cur = dbhandle.cursor()
    # Execute
    cur.execute(query)
    tuples = cur.fetchall()

    # Return rec array
    names = [d[0] for d in cur.description]

    if not tuples:
        types = []
        for d in cur.description:
            types.append((d[0], object))
        return numpy.recarray((0,), dtype=types)
    return numpy.rec.array(tuples, names=names)

###################################################################
def get_coadd_from_pfwid(dbi, pfwid, band):
    """ Method to fing the coadd file name from the pfw attempt id

        Parameters
        ----------
        dbi - database handle

        pfwid - int
            The pfw attempt id

        band - str
            The current band

        Returns
        -------
        The coadd file name associated with the pfw attempt id and band

    """
    curs = dbi.cursor()
    curs.execute(f"select filename from desfile where pfw_attempt_id={int(pfwid):d} and filetype='coadd' and (filename like '%_{band}.%' or filename like '%_{band}_%') and compression is null")
    try:
        results = curs.fetchone()[0]
    except:
        print(f"Could not find coadd file for pfw_attempt_id: {int(pfwid):d}  band: {band}")
        raise
    return results

###################################################################
def get_pfwid_from_coadd(dbi, coadd_file):
    """ Method to find the pfw attempt id from the coadd file name

        Parameters
        ----------
        dbi - database handle

        coadd_file - str
            The name of the coadd file to use

        Returns
        -------
        The pfw attempt id the generated the coadd file
    """
    curs = dbi.cursor()
    tcoadd_file = coadd_file.split('/')[-1]
    if tcoadd_file.endswith('.fits'):
        temp = f"'{tcoadd_file}' and compression is null"
    else:
        parts = tcoadd_file.split('.fits')
        temp = f"'{parts[0]}.fits' and compression='{parts[1]}'"
    sql = f"select pfw_attempt_id from desfile where filename={temp}"

    curs.execute(sql)
    results = curs.fetchall()
    if len(results) != 1:
        raise Exception("Could not determine the pfw_attempt_id from the coadd_file")
    return results[0][0]

###################################################################
def get_nwgint_files(dbi, band, pfwid=None, coadd_file=None):
    """ Method to get nwgint file names by querying the DB

        Parameters
        ----------
        dbi - database handle

        band - str
            The band being processed

        pfwid - int
            The pdw_attampt_id being processed

        Returns
        -------
        Tuple containing a dictionary of the initial header data for the nwgint files, and a list of the file names

    """
    # query for the file names
    curs = dbi.cursor()
    if coadd_file is not None:
        tcoadd_file = coadd_file.split('/')[-1]
        if tcoadd_file.endswith('.fits'):
            temp = f"'{tcoadd_file}' and compression is null"
        else:
            parts = tcoadd_file.split('.fits')
            temp = f"'{parts[0]}.fits' and compression='{parts[1]}'"
        sql = f"select filename from desfile where id in(select parent_desfile_id from opm_was_derived_from start with child_desfile_id in (select id from desfile where filename={temp}) connect by nocycle prior parent_desfile_id=child_desfile_id and level<5)and filetype='coadd_nwgint'"
    else:
        sql = f"select filename from image where pfw_attempt_id={int(pfwid):d} and band='{band}' and filetype='coadd_nwgint'"

    curs.execute(sql)
    results = curs.fetchall()
    listinfo = {}
    filenames = []

    for res in results:
        filenames.append(res[0])
        listinfo[res[0]] = {'fullname': res[0]}

    return listinfo, filenames

###################################################################
def get_mag_zero(dbi, filename, table, zsource, zversion, z2source, z2version, zpt2, gray_tbl):
    """ get mag_zero and sigma_mag_zero"""

    curs = dbi.cursor()
    curs.execute(f"select MAG_ZERO, SIGMA_MAG_ZERO, EXPNUM, CCDNUM from {table} where imagename='{filename}' and source='{zsource}' and version='{zversion}'")
    res = curs.fetchone()
    if res is not None:
        mag_zero, sigma_mag_zero, expnum, ccd = res[0], res[1], res[2], res[3]
    elif z2source is not None:
        curs.execute(f"select MAG_ZERO, SIGMA_MAG_ZERO, EXPNUM, CCDNUM from {zpt2} where imagename='{filename}' and source='{z2source}' and version='{z2version}'")
        res = curs.fetchone()
        if res is None:
            raise Exception(f"Cannot locate a mag_zero value for {filename}")
        mag_zero, sigma_mag_zero, expnum, ccd = res[0], res[1], res[2], res[3]
    else:
        raise Exception(f"Cannot locate a mag_zero value for {filename}")

    curs.execute(f"select FGCM_GRY from {gray_tbl} where expnum={int(expnum):d} and ccdnum={int(ccd):d} and pfw_attempt_id in (select pfw_attempt_id from desfile where filename='{filename}' and compression is null)")
    res = curs.fetchone()
    if res is not None:
        return mag_zero, sigma_mag_zero, res[0]
    return mag_zero, sigma_mag_zero, 0.

###################################################################
def get_redimg_files(dbi, nwgtab, table, zsource, zversion, zflag, im_to_tile, tilename, coadd_file=None, z2source=None,
                     z2version=None, z2flag=None, zpt2=None, gry_tbl=None):
    """ Method to retrieve that red image header data from the database

        Parameters
        ----------
        dbi - database handle

        nwgtab - dict
            Dict containing the data of the nwgint header data

        table - str
            The name of the zaero point table to use

        zsource - str
            The source to use for the zero point data (matches entry in source column of zero point table)

        zversion - str
            The version of the z source to use (matches entry in the version column of the zero point table)

        im_to_tile - str
            The table to use to locate the correct red images for the specified tile

        tilename - str
            The name of the tile currently being processed.

        Returns
        -------
        Tuple containing a dictionary of the red img header data and an list of the file names
    """
    listinfo = {}
    filenames = []
    curs = dbi.cursor()
    filelist = []

    if coadd_file is not None:
        tcoadd_file = coadd_file.split('/')[-1]
        if tcoadd_file.endswith('.fits'):
            temp = f"'{tcoadd_file}' and compression is null"
        else:
            parts = tcoadd_file.split('.fits')
            temp = f"'{parts[0]}.fits' and compression='{parts[1]}'"
        sql = f"select distinct(filename) from desfile where id in(select parent_desfile_id from opm_was_derived_from start with child_desfile_id in (select id from desfile where filename={temp}) connect by nocycle prior parent_desfile_id=child_desfile_id and level<8)and filetype='red_immask'"
        curs.execute(sql)
        results = curs.fetchall()
        for res in results:
            filelist.append(res[0])
    else:
        # for each of the nwgint file names get the corresponding red image header info
        for i in range(len(nwgtab['FILENAME'])):
            sql = f"select filename from image where expnum={nwgtab['EXPNUM'][i]:d} and ccdnum={nwgtab['CCDNUM'][i]:d} and filetype='red_immask' and filename in (select filename from {im_to_tile} where tilename='{tilename}')"

            curs.execute(sql)
            results = curs.fetchall()
            if len(results) != 1:
                raise Exception(f"Incorrect number of files {len(results)}")

            filelist.append(results[0][0])
        # get the mag_zero value
    for filename in filelist:
        #curs.execute("select MAG_ZERO from %s where imagename='%s' and source='%s' and version='%s'" %
        #             (table, filename, zsource, zversion))
        #mag_zero = curs.fetchone()
        #if mag_zero is not None:
        #    mag_zero = mag_zero[0]
        #elif z2source is not None:
        #    curs.execute("select MAG_ZERO from %s where imagename='%s' and source='%s' and version='%s'" %
        #                 (zpt2, filename, z2source, z2version))
        #    mag_zero = curs.fetchone()
        #    if mag_zero is None:
        #        raise Exception("Cannot locate a mag_zero value for %s" % filename)
        #    else:
        #        mag_zero = mag_zero[0]
        #else:
        #    raise Exception("Cannot locate a mag_zero value for %s" % filename)
        mag_zero, sigma_mag_zero, fgcm_gry = get_mag_zero(dbi, filename, table, zsource, zversion, z2source, z2version, zpt2, gry_tbl)
        filenames.append(filename)
        listinfo[filename] = {'fullname': filename,
                              'mag_zero': float(mag_zero),
                              'sigma_mag_zero': float(sigma_mag_zero),
                              'fgcm_gry': float(fgcm_gry)}
    return listinfo, filenames

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

    sql = """select i.*, pq.fwhm_mean*0.263 as FWHM, pq.FWHM_FROMFLUXRADIUS_MEAN*0.263 as FWHM_FLUXRAD,
             case when qa.t_eff < 0.0 then qa.f_eff*qa.b_eff
             else qa.t_eff end as T_EFF from image i, opm_filename_gtt g,miscfile m,
             psf_qa pq, qa_summary qa where i.filename=g.filename and m.pfw_attempt_id=i.pfw_attempt_id and
             m.ccdnum=i.ccdnum and m.filename=pq.filename and m.filetype='xml_psfex' and qa.expnum=i.expnum
             order by i.filename"""


    if miscutils.fwdebug_check(3, "MANGLEDB_DEBUG"):
        miscutils.fwdebug_print(f"sql = {sql}")

    red_tab = despyastro.query2dict_of_columns(sql, dbi, True)
    red_tab["T_EFF_EXPTIME"] = red_tab["T_EFF"] * red_tab["EXPTIME"]
    return red_tab


###################################################################
def get_streak_info(dbi, schema='', redfiles=None):
    """ Read streak information from database
    """

    if redfiles is not None:
        load_gtt_filename(dbi, redfiles)

    sql = """select i.filename as redfilename, s.*
             from streak s, catalog c, image i, opm_filename_gtt g where
             s.filename=c.filename and i.filename=g.filename and
             i.expnum=c.expnum and i.ccdnum=c.ccdnum and i.pfw_attempt_id=c.pfw_attempt_id and not
             (s.ra_1=s.ra_2 and s.ra_1=s.ra_3 and s.ra_1=s.ra_4 and s.dec_1=s.dec_2 and s.dec_1=s.dec_3 and s.dec_1=s.dec_4)
          """
#   RAG 20190204: added constraint to remove corner case streaks that have zero size
#       (i.e. ra1=ra2=ra3=ra4 and dec1=dec2=dec3=dec4)
#   This corner case appears to have arisen from the addition of the streak connector in single-epoch processing
#   (speculating... for streaks with grazing incidence on another chip)

    if miscutils.fwdebug_check(3, "MANGLEDB_DEBUG"):
        miscutils.fwdebug_print(f"sql = {sql}")

    streak_tab = query2rec(sql, dbhandle=dbi)
    return streak_tab


###################################################################
def get_satstar_info(dbi, schema='', redfiles=None, max_radius=400.0):
    """ Read satstar information from database """

    if redfiles is not None:
        load_gtt_filename(dbi, redfiles)

    sql = f"""select i.filename as redfilename, s.*
             from satstar s, catalog c, image i, opm_filename_gtt g where
             s.filename=c.filename and s.radius<{max_radius:f} and s.radius>0.1 and i.filename=g.filename and
             i.expnum=c.expnum and i.ccdnum=c.ccdnum and i.pfw_attempt_id=c.pfw_attempt_id
          """

    if miscutils.fwdebug_check(3, "MANGLEDB_DEBUG"):
        miscutils.fwdebug_print(f"sql = {sql}")

    satstar_tab = query2rec(sql, dbhandle=dbi)
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
        miscutils.fwdebug_print(f"sql = {sql}")

    bleedtrail_tab = query2rec(sql, dbhandle=dbi)
    return bleedtrail_tab

###################################################################
def get_headfile_info(dbi, me_attid, redfiles=None):
    """ Read coadd head file information from database """

    if redfiles is not None:
        load_gtt_filename(dbi, redfiles)

    sql = f"""select i.filename as redfilename, m.filename
             from miscfile m, image i, opm_filename_gtt g where
             i.filename=g.filename and
             i.expnum=m.expnum and i.ccdnum=m.ccdnum and m.pfw_attempt_id={me_attid}
             and m.filetype='coadd_head_scamp'
          """

    if miscutils.fwdebug_check(3, "MANGLEDB_DEBUG"):
        miscutils.fwdebug_print(f"sql = {sql}")

    headfile_tab = despyastro.query2dict_of_columns(sql, dbi, True)
    return headfile_tab

##################################################################
def get_nwgint_info_from_reds(dbi, me_attid, redfiles=None):
    """ Read nwgint information from database """

    if redfiles is not None:
        load_gtt_filename(dbi, redfiles)

    sql = f"""select i1.filename as redfilename, i2.*
             from image i1, image i2, opm_filename_gtt g where
             i1.filename=g.filename and
             i1.expnum=i2.expnum and i1.ccdnum=i2.ccdnum and i2.pfw_attempt_id={me_attid}
             and i2.filetype='coadd_nwgint'
          """

    if miscutils.fwdebug_check(3, "MANGLEDB_DEBUG"):
        miscutils.fwdebug_print(f"sql = {sql}")

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
        miscutils.fwdebug_print(f"sql = {sql}")

    nwgint_tab = despyastro.query2dict_of_columns(sql, dbi, True)
    return nwgint_tab


###################################################################
def get_tileid(tilename, db_section=None, des_services=None):
    """ Get tile id from database for given tilename """

    dbh = despydb.desdbi.DesDbi(des_services, db_section, True)

    # Prepare the query
    query = f"select id from coaddtile_geom where tilename='{tilename}'"

    # Get the cursor
    cur = dbh.cursor()
    cur.execute(query)

    item = cur.fetchone()
    return item[0]


###################################################################
def get_coadd_object_info(dbi, tilename, pfw_attempt_id,
                          ra_column, dec_column, dbschema=''):
    """ Get id and location information for coadd objects """


    sql = f"""select co.id, co.filename, co.object_number as "NUMBER",
                    co.{ra_column}, co.{dec_column}
             from {dbschema}coadd_object co, {dbschema}catalog c
             where co.filename=c.filename and
                   c.filetype='coadd_det_cat' and
                   c.tilename='{tilename}' and
                   c.pfw_attempt_id={pfw_attempt_id}
          """

    if miscutils.fwdebug_check(3, "MANGLEDB_DEBUG"):
        miscutils.fwdebug_print(f"sql = {sql}")

    coadd_object_tab = query2rec(sql, dbhandle=dbi)
    return coadd_object_tab

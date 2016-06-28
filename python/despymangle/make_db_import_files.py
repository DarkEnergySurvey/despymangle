#!/usr/bin/env python

# $Id$
# $Rev::                                  $:  # Revision of last commit.
# $LastChangedBy::                        $:  # Author of last commit.
# $LastChangedDate::                      $:  # Date of last commit.

import despymangle.mangle_db as mdb

def make_csv_files(config, dbi):

    coadd_object_tab = mdb.get_coadd_object_info(dbi, 
                                                 config['tilename'],
                                                 config['pfw_attempt_id'],
                                                 config['schema'], 
                                                 config['ra_column'],
                                                 config['dec_column'])

    #print coadd_object_tab
    #print coadd_object_tab['RA']
    #print coadd_object_tab['DEC']
    #print coadd_object_tab['ID']


    fndb_pattern = config['outputdir'] + config['dbcsv'] + '_%s.csv'
    fndb_ccdgon = fndb_pattern % ('ccdgon')
    fndb_ccdmoly = fndb_pattern % ('ccdmoly')
    fndb_molygon_ccdgon = fndb_pattern % ('molygon_ccdgon')
    fndb_coadd_object_molygon = fndb_pattern % ('coadd_object_molygon') 

    print 'fndb_ccdgon =', fndb_ccdgon
    print 'fndb_ccdmoly =', fndb_ccdmoly
    print 'fndb_molygon_ccdgon =', fndb_molygon_ccdgon
    print 'fndb_coadd_object_molygon =', fndb_coadd_object_molygon



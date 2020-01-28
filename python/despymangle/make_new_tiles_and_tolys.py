#!/usr/bin/env python3

import os
import sys
import numpy as np
import set_basic_mangle_params_Y3A1v1 as par

import despydb
import despymangle.mangle_utils as mu
import pymangle as pym



# Get DB hadle:
# The try/except is only needed in case your credentials live in a
# file different than: $HOME/.desservices.ini

try:
    desdmfile = os.environ["des_services"]
except KeyError:
    desdmfile = None

section = "db-desoper"
dbh = despydb.desdbi.DesDbi(desdmfile,section)



### get the tiles file

# Prepare the query


##query="""select coaddtile_id, raul, decul, rall, decll,  ralr, declr, raur, decur from coaddtile  """


###query = """ select c.*, cc.coaddtile_id from felipe.coaddtile_new c, coaddtile cc where c.tilename= cc.tilename """
query = "select c.id, c.rac1, c.decc1, c.rac2, c.decc2, c.rac3, c.decc3, c.rac4, c.decc4, c.uramin, c.uramax, c.udecmin, c.udecmax, c.crossra0 from felipe.coaddtile_geom c"




# Get the cursor
cur = dbh.cursor()
cur.execute(query)


items = cur.fetchall()
#itemprint items


AA = np.array(items, dtype=np.dtype('a25'))
g = open('jtiles', 'w')
N = len(items)

#for i in range(len(items)):
for i in range(N):
#    print i
    id = AA[i, 0]
    ll = np.float128(AA[i, 1: -1])
    # if ( (ll[1]<-30)*(ll[1]>-32)):
    print(i)
    print(f"edges  {%id}", file=g)
    aa = mu.c2e(ll[0], ll[1], ll[6], ll[7], ll[4], ll[5],ll[2], ll[3])
    print(aa[1], file=g)

g.close()



if 1 == 1:
    os.system(f"poly2poly  -ie2 {par.mtol} jtiles 'test_tiles.pol")
    os.system(f"snap -S {par.snap} {par.mtol} test_tiles.pol stest_tiles.pol")
    os.system(f"pixelize {par.pix} {par.mtol} stest_tiles.pol pstest_tiles.pol")
    os.system(f"snap {par.snap} {par.mtol} pstest_tiles.pol spstest_tiles.pol")
    os.system(f"balkanize -Ba   {par.mtol} spstest_tiles.pol bspstest_tiles.pol")

cur.close()



## get tolys

g = open('jtolys', 'w')
for i in range(N):
    print(i)
    id = AA[i, 0]
    ll = np.float128(AA[i, 1: -1])
    crossra = AA[i, -1]
    #if ( (ll[1]<-30)*(ll[1]>-32)):
    print(i)
    print(f"edges  {id}", file=g)
    if crossra == 'Y':
        print('youpi')
        aa = mu.c2e(ll[9], ll[10], ll[8], ll[10], ll[8], ll[11], ll[9], ll[11])
    else:
        aa = mu.c2e(ll[8], ll[10], ll[9], ll[10], ll[9], ll[11], ll[8], ll[11])
    print(aa[1], file=g)

g.close()
if 1 == 1:
    os.system(f"poly2poly  -ie2 {par.mtol} jtolys test_tolys.pol")
    os.system(f"snap -S {par.snap} {par.mtol} test_tolys.pol stest_tolys.pol")
    os.system(f"pixelize {par.pix} {par.mtol} stest_tolys.pol pstest_tolys.pol")
    os.system(f"snap {par.snap} {par.mtol} pstest_tolys.pol spstest_tolys.pol")
    os.system(f"balkanize -Ba -vo {par.mtol} spstest_tolys.pol bspstest_tolys.pol")

#### make tiles

if 1 == 1:
    new_tile_dir = par.tiledir
    if not os.path.isdir(new_tile_dir):
        os.mkdir(new_tile_dir)
    ##'/home/benoitl/DESdata/tiles_test/DES/OPS/'
    #tiledir=

    os.system(f"cp stest_tiles.pol {new_tile_dir}{par.project}_tiles_unpix.pol")
    os.system(f"cp stest_tolys.pol {new_tile_dir}{par.project}_tolys_unpix.pol")

    os.system(f"cp spstest_tiles.pol {new_tile_dir}{par.project}_tiles{par.restag}.pol")
    os.system(f"cp spstest_tolys.pol {new_tile_dir}{par.project}_tolys{par.restag}.pol")

    os.system(f"poly2poly -od {par.mtol} spstest_tiles.pol {par.project}_tiles{par.restag}.dpol")
    os.system(f"poly2poly -od {par.mtol} spstest_tolys.pol {par.project}_tolys{par.restag}.dpol")

    os.system(f"mv {par.project}_tiles{par.restag}.dpol {new_tile_dir}")
    os.system(f"mv {par.project}_tolys{par.restag}.dpol {new_tile_dir}")


    if 1 == 1:
        os.system('rm jtiles test_tiles.pol  stest_tiles.pol  pstest_tiles.pol  spstest_tiles.pol bspstest_tiles.pol')
        os.system('rm jtolys test_tolys.pol  stest_tolys.pol  pstest_tolys.pol  spstest_tolys.pol bspstest_tolys.pol')

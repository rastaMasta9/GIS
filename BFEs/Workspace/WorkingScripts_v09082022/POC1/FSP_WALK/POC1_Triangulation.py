import time
s = time.time()
from new_methods import *
import geopandas as gpd
import pandas as pd


bfe = gpd.read_file('bfe_mre_all.shp')
fsp = gpd.read_file('fsp_sample_mre1.shp')

# Reproject files to UTM
bfe.to_crs(26913, inplace=True)
bfe = bfe[['ELEV', 'geometry']]

fsp.to_crs(26913, inplace=True)
fsp = fsp[['geometry']]
fsp.reset_index(inplace=True)

# Correct broken BFEs if there is any
bfe_bfe = bfe.sjoin(bfe, how='left')
bfe_bfe.reset_index(inplace=True)
brok = bfe_bfe.loc[bfe_bfe['index'] != bfe_bfe['index_right']]
if brok.shape[0] > 0:
    brok_geoms = brok[['index', 'geometry']]
    fix = brok.merge(brok_geoms, left_on='index_right', right_on='index')
    fix['new_geom'] = fix.apply(lambda x: linemerge(list((x.geometry_x, x.geometry_y))), axis=1)
    fix_gpd = gpd.GeoDataFrame(fix[['index_x', 'ELEV_left', 'new_geom']], geometry='new_geom', crs=26913)
    fix_gpd.rename(columns={'index_x': 'index', 'ELEV_left': 'ELEV', 'new_geom': 'geometry'}, inplace=True)
    fix_gpd.set_index('index', inplace=True)
    bfe_clean = bfe.loc[~bfe.index.isin(fix_gpd.index)]
    bfe = pd.concat([bfe_clean, fix_gpd])
    bfe = bfe.drop_duplicates('geometry')
   

# Extend BFEs over FSP Poly then reset bfe to new position
bfe_extend = extend_bfe(bfe, fsp)
bfe = bfe_extend.clip(fsp)

# Split FSP Poly by extended BFEs
fsp_s = split_fsp(fsp, bfe_extend)

triangles = gpd.GeoDataFrame()
# Iterate through each FSP polygon
for f in fsp_s['geometry']:
    f = g(f, 26913)
    bfe_set = bfe.sjoin(f)
    bfe_set = bfe_set[['ELEV', 'geometry']]
    bfe_set = remove_multiline_BFE(bfe_set)
    bfe_count = bfe_set.shape[0]
    
    # Getting Z-geom for BFE Points
    bfe_pts = bfe_zpts(bfe_set)[['ELEV', 'geometry']]
    
    
    # FSP Simplify, Interpolation, and Z-geom
    fsp_i_pts = (f.pipe(fsp_pts_simplify, tolerance=3)
                .pipe(IDW, bfe=bfe_set, power=2, bfe_count=bfe_count)
                .pipe(ELEV_2geom)
    )

    # Concat and Triangulate
    all_pts = pd.concat([bfe_pts, fsp_i_pts], ignore_index=True)
    all_pts_multigeom = MultiPoint(all_pts.geometry.to_list())

    tin = triangulate(all_pts_multigeom)
    tin_df = g(tin, 26913)

    # Extract Geom
    final_tin = extract_geom(tin_df)
    triangles = pd.concat([triangles, final_tin], ignore_index=True)

triangles.overlay(fsp).to_file(input('Output File: '))
eta = time.time()
print('Processing Complete: ', (eta - s))
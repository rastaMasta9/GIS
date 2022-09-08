from methods import *
import geopandas as gpd
import pandas as pd

bfe = gpd.read_file('bfe_mre_all.shp')
fsp = gpd.read_file('fsp_sample_mre1.shp')
fl = gpd.read_file('Flowline_sample.shp')

# Reproject files to UTM
bfe.to_crs(26913, inplace=True)
bfe = bfe[['ELEV', 'geometry']]

fsp.to_crs(26913, inplace=True)
fsp = fsp[['geometry']]
fsp.reset_index(inplace=True)

fl.to_crs(26913, inplace=True)
fl = fl[['geometry']]
fl.reset_index(inplace=True)

# Extend BFEs over FSP Poly
bfe_extend = extend_bfe(bfe, fsp)

# Split FSP Poly by extended BFEs
fsp_s = split_fsp(fsp, bfe_extend)

# Split Flowline by BFEs
fl = split_flowline(fl, bfe)

# Create Flowline Buffer for BFE intersection
# 'fl_buff' = Original Flowline buffer
# 'fl_buff_sj' = Flowline Buffer with Adjacent Flowlines
fl_buff, fl_buff_sj = buff_fl_intersection(fl, 1)

# 'adj' groupby adjacent flowline seg.IDs
adj = create_adjacent_gpd(fl_buff_sj)

# Filter 'adj' for Forks( i.e. FLowlines with more than 2 intersections)
#adj['threes'] = adj.apply(lambda x: 1 if len(x['index_right']) == 3 else np.nan, axis=1)
#threes = adj.loc[adj['threes'].notnull()]

# List of groups of Forks
#if threes.shape[0] > 0:
#    forks_list = find_forks(threes)

# Run Triangulation Build by Walking FLowlines with 2 BFE intersections
triangles = gpd.GeoDataFrame()
for i, a in adj.iterrows():
    if len(a['index_right']) == 1: 
        pass
    elif len(a['index_right']) == 3:
        pass
    else:
        fl_i = fl.iloc[i]
        fl_buff_i = fl_buff.loc[fl_buff['index'] == i]

        bfe_set = bfe.sjoin(fl_buff_i, how='left', predicate='intersects')
        bfe_set = bfe_set.loc[bfe_set['index_right'].notnull()]
        bfe_set = bfe_set[['index', 'ELEV', 'geometry']]

        # Getting Z-geom for BFE Points
        bfe_pts = bfe_zpts(bfe_set)

        # BFE Centroid and Flowline interpolation
        fl_i_pts = flowline_interpolation(bfe_set, fl_i, 15)

        # FSP Simplify
        fsp_pts = fsp_pts_simplify(fsp_s, fl_i_pts, tolerance=3)

        # FSP Interpolation
        fsp_i_pts = IDW(bfe_set, fsp_pts)

        # Concat and Triangulate
        all_pts = pd.concat([fl_i_pts, fsp_i_pts], ignore_index=True)
        all_pts_multigeom = MultiPoint(all_pts.geometry.to_list())

        tin = triangulate(all_pts_multigeom)
        tin_df = g(tin, 26913)
        
        # Extract Geom
        final_tin = extract_geom(tin_df)
        triangles = pd.concat([triangles, final_tin], ignore_index=True)

triangles.to_file(input('Output File: '))






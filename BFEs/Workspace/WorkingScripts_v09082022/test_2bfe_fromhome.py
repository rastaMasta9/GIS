from methods import *
import geopandas as gpd
import pandas as pd

bfe = gpd.read_file('bfe_mre_all.shp')
fsp = gpd.read_file('fsp_sample_mre1.shp')
fl = gpd.read_file('Flowline_Forks.shp')

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
fl['buff'] = fl.apply(lambda x: x.geometry.buffer(1))
fl_buff = fl.set_geometry('buff')
fl_buff = fl_buff[['buff']]
fl = fl[['geometry']]

# Find Forks
forks_ls = find_forks(fl)

# Parse nonforks df out of main
all_fsegs_forks = [segs for e in forks_ls for segs in e]
fl_non_forks = fl.loc[~fl.index.isin(all_fsegs_forks)].reset_index()
fl_buff_non_forks = fl_buff.loc[~fl_buff.index.isin(all_fsegs_forks)].reset_index()

# Run Triangulation Build by Walking FLowlines with 2 BFE intersections
triangles = gpd.GeoDataFrame()
fingers = []
for i, f in fl_non_forks.iterrows():
    if len(a['index_right']) != 2: 
        pass
    else:
        fl_i = fl.iloc[i]
        fl_buff_i = fl_buff.loc[fl_buff['index'] == i]

        bfe_set = bfe.sjoin(fl_buff_i, how='left', predicate='intersects')
        bfe_set = bfe_set.loc[bfe_set['index_right'].notnull()]
        bfe_set = bfe_set[['ELEV', 'geometry']]

        if bfe_set.shape[0] != 0: 
            # Getting Z-geom for BFE Points
            bfe_pts = bfe_zpts(bfe_set)

            # BFE Centroid and Flowline interpolation
            fl_i_pts = flowline_interpolation(bfe_set, fl_i, divisions=15, power=2)

            # FSP Simplify
            fsp_pts = fsp_pts_simplify(fsp_s, fl_i_pts, tolerance=3)

            # FSP Interpolation
            fsp_i_pts = IDW(bfe_set, fsp_pts, power=2)

            # Concat and Triangulate
            all_pts = pd.concat([bfe_pts, fl_i_pts, fsp_i_pts], ignore_index=True)
            all_pts_multigeom = MultiPoint(all_pts.geometry.to_list())

            tin = triangulate(all_pts_multigeom)
            tin_df = g(tin, 26913)
            
            # Extract Geom
            final_tin = extract_geom(tin_df)
            triangles = pd.concat([triangles, final_tin], ignore_index=True)
            
        else:
            fingers.append(i)
#fl = fl.loc[~fl.index.isin(fingers)]
#fl_buff = fl_buff.loc[~fl_buff.index.isin(fingers)]
        
# fork workflow
if forks_list:
    fork_triangles = gpd.GeoDataFrame()
    for f in forks_list:
        
        f_seg, buff_union = union_fork(f, fl, fl_buff)
        bfe_set = bfe.sjoin(buff_union, how='left', predicate='intersects')
        bfe_set = bfe_set.loc[bfe_set['index_right'].notnull()]
        bfe_set = bfe_set[['ELEV', 'geometry']]

        # Getting Z-geom for BFE Points
        bfe_pts = bfe_zpts(bfe_set)

        # Flowline interpolation
        if bfe_set.shape[0] == 2:
            fake_fork = True
        else:
            fake_fork = False
            
        f_interp_df = flowline_interpolation(bfe_set, f_seg, divisions=15, power=2, fork='yes', non_fork=fake_fork)
            
        print('bfe_set')
        print(bfe_set)
        print('___________')
        print('fork')
        print(f)
        print('____________')
        print(f_interp_df)
        # FSP Simplify and interpolation
        fsp_pts = fsp_pts_simplify(fsp_s, f_interp_df, tolerance=3)
        fsp_interp_df = IDW_Forks(bfe_set, fsp_pts, power=2, non_fork=fake_fork)

        # Concat and Triangulate
        all_pts = pd.concat([bfe_pts, f_interp_df, fsp_interp_df], ignore_index=True)
        all_pts_multigeom = MultiPoint(all_pts.geometry.to_list())

        tin = triangulate(all_pts_multigeom)
        tin_df = g(tin, 26913)
        
        # Extract Geom
        final_tin = extract_geom(tin_df)
        fork_triangles = pd.concat([fork_triangles, final_tin], ignore_index=True)


triangles.to_file('triangles09_4.shp') 
fork_triangles.to_file('fork_triangles09_4.shp')    
"""
print('Triangles')
print(triangles['geometry'].value_counts())
print('_____________')
print('Forks')
print(fork_triangles['geometry'].value_counts())
"""
#all_triangles = pd.concat([triangles, fork_triangles], ignore_index=True)

#all_triangles.to_file(input('Output File: '))

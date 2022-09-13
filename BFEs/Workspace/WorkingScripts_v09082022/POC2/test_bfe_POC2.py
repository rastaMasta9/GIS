import time
s = time.time()
from methods_fromhome import *
import geopandas as gpd
import pandas as pd



bfe = gpd.read_file('bfe_mre_all_poc2.shp')
fsp = gpd.read_file('fsp_2.shp')
fl = gpd.read_file('flowline_poc2.shp')

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

# Find forks and place in list
forks_ls = find_forks(fl)

# Create Flowline Buffer
fl['buff'] = fl.apply(lambda x: x.geometry.buffer(1))
fl_buff = fl.set_geometry('buff')
fl_buff = fl_buff[['buff']]
fl = fl[['geometry']]

# get fl_segs that are forks
all_fsegs_forks = [segs for e in forks_ls for segs in e]

fl_non_forks = fl.loc[~fl.index.isin(all_fsegs_forks)].reset_index()
fl_buff_non_forks = fl_buff.loc[~fl_buff.index.isin(all_fsegs_forks)].reset_index()

triangles = gpd.GeoDataFrame()
for i, f in fl_non_forks.iterrows():
    fl_buff_i = fl_buff_non_forks.loc[fl_buff_non_forks['index'] == f['index']]
    
    bfe_set = bfe.sjoin(fl_buff_i, how='left', predicate='intersects')
    bfe_set = bfe_set.loc[bfe_set['index_right'].notnull()]
    bfe_set = bfe_set[['ELEV', 'geometry']]
    if bfe_set.shape[0] == 1:
        pass
    else:
        # Getting Z-geom for BFE Points
        bfe_pts = bfe_zpts(bfe_set)
        print(bfe_set)
        # BFE Centroid and Flowline interpolation
        fl_i_pts = flowline_interpolation(bfe_set, 
                                            f,
                                            line_segs=None,
                                            buff_segments=None, 
                                            divisions=15, 
                                            power=2
        )
                                            
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

fork_triangles = gpd.GeoDataFrame()
for f in forks_ls:
    f_segs, f_union, buff_segs, buff_union = union_fork(f, fl, fl_buff)
    bfe_set = bfe.sjoin(buff_union, how='left', predicate='intersects')
    bfe_set = bfe_set.loc[bfe_set['index_right'].notnull()]
    bfe_set = bfe_set[['ELEV', 'geometry']]
    
    # Getting Z-geom for BFE Points
    bfe_pts = bfe_zpts(bfe_set)

    if bfe_set.shape[0] == 2:
        fake_fork = True
    else:
        fake_fork = False

    f_interp_df = flowline_interpolation(bfe_set, 
                                        f_union,
                                        f_segs,
                                        buff_segs, 
                                        divisions=15, 
                                        power=2, 
                                        fork='yes', 
                                        non_fork=fake_fork
    )

    # FSP Simplify and interpolation
    fsp_pts = fsp_pts_simplify(fsp_s, f_interp_df, tolerance=3)
    fsp_interp_df = IDW_Forks(bfe_set, 
                                fsp_pts, 
                                power=2, 
                                non_fork=fake_fork
    )
    
    # Concat and Triangulate
    all_pts_f = pd.concat([bfe_pts, f_interp_df, fsp_interp_df], ignore_index=True)

    all_pts_multigeom_f = MultiPoint(all_pts_f.geometry.to_list())

    tin_f = triangulate(all_pts_multigeom_f)
    tin_df_f = g(tin_f, 26913)
    
    # Extract Geom
    final_tin_f = extract_geom(tin_df_f)
    fork_triangles = pd.concat([fork_triangles, final_tin_f], ignore_index=True)

all_triangles = pd.concat([triangles, fork_triangles], ignore_index=True)
end = (time.time() - s)
all_triangles.to_file(input('Enter output file: '))
print(f'Processing Time: {end}')




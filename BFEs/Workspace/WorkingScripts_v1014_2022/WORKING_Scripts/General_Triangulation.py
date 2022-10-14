import time
s = time.time()
from new_methods import *

bfe = gpd.read_file(input('Enter BFE shapefile: '))
fsp = gpd.read_file(input('Enter Flood Zone Shapefile: '))
out = input('Enter Prefix for all output files: ') # County FIPS
# Reproject files to UTM
bfe.to_crs(26913, inplace=True)
bfe = bfe[['ELEV', 'geometry']]

fsp.to_crs(26913, inplace=True)
fsp = fsp[['geometry']]
fsp.reset_index(inplace=True)

# Correct broken BFEs if there is any
bfe_bfe = bfe.sjoin(bfe, how='left')
bfe_bfe.reset_index(inplace=True)

# Indices will not match if seperate LineStrings join
brok = bfe_bfe.loc[bfe_bfe['index'] != bfe_bfe['index_right']]
if brok.shape[0] > 0:
    brok_geoms = brok[['index', 'geometry']]

    # Merge df with correct geometries then merge together
    fix = brok.merge(brok_geoms, left_on='index_right', right_on='index')
    fix['new_geom'] = fix.apply(lambda x: linemerge(list((x.geometry_x, x.geometry_y))), axis=1)

    # CLean, rename df
    fix_gpd = gpd.GeoDataFrame(fix[['index_x', 'ELEV_left', 'new_geom']], geometry='new_geom', crs=26913)
    fix_gpd.rename(columns={'index_x': 'index', 'ELEV_left': 'ELEV', 'new_geom': 'geometry'}, inplace=True)
    fix_gpd.set_index('index', inplace=True)

    # Remove orginal BFEs and replace with new merged BFEs
    bfe_clean = bfe.loc[~bfe.index.isin(fix_gpd.index)]
    bfe = pd.concat([bfe_clean, fix_gpd])
    bfe = bfe.drop_duplicates('geometry')


# Extend BFEs over FSP Poly then reset bfe to new position
print('Extending BFE...')
bfe_ext = bfe_extend(bfe, fsp)
print('Clipping BFE...')
bfe = bfe_ext.overlay(fsp)

print('Splitting FSP...')
# Split FSP Poly by extended BFEs
fsp_s = split_fsp(fsp, bfe_ext)

# Remove slivers due to overlaps from extended BFEs
fsp_s = fsp_s.overlay(fsp)
fsp_s = fsp_s.loc[fsp_s.geom_type == 'Polygon']
fsp_s = fsp_s[['geometry']]
fsp_s.reset_index(inplace=True)

print('Splitting FSP...')
# Split FSP Poly by extended BFEs
fsp_s = split_fsp(fsp, bfe_ext)

# Remove slivers due to overlaps from extended BFEs
print('Clipping FSP extent back to original...')
fsp_s = fsp_s.overlay(fsp)
fsp_s = fsp_s.loc[fsp_s.geom_type == 'Polygon']
fsp_s = fsp_s[['geometry']]
fsp_s.reset_index(inplace=True)

# Iterate through each FSP polygon
triangles = gpd.GeoDataFrame()
poly_errors = gpd.GeoSeries(crs=4326)
for i, f in fsp_s.iterrows():
    f = g(f.geometry, 26913)
    bfe_set = bfe.sjoin(f)
    
    bfe_set = bfe_set[['ELEV', 'geometry']]
    bfe_set = remove_multiline_BFE(bfe_set)
    
    print('POLY INDEX: ', i)
    #Ignoring Potential Polygon Slivers due to extension of BFEs overlapping small poritions of FSP
    if len(f.sjoin(bfe_set)) == 0:
        continue
    # Potential Island poly which will not intersect BFEs. Ignore!
    elif bfe_set.shape[0] == 0:
        continue
    else:
        
        # Getting Z-geom for BFE Points
        bfe_pts_utm, bfe_pts_84 = bfe_zpts(bfe_set)

        # FSP Simplify, Interpolation, and Z-geom
        fsp_i_pts = fsp_pts(f, bfe_pts=bfe_pts_utm, bfe_set=bfe_set, diff_area=800)
        
        if fsp_i_pts is not None:
            print('Creating more points..')
            fsp_i_pts = (fsp_i_pts.pipe(IDW, bfe=bfe_set, power=2)
                        .pipe(ELEV_2geom)
        )

            # Concat and Triangulate
            all_pts = pd.concat([bfe_pts_84, fsp_i_pts], ignore_index=True)
            all_pts_multigeom = MultiPoint(all_pts.geometry.to_list())

            tin = triangulate(all_pts_multigeom)
            tin_df = g(tin, 4326)
        
            # Extract Geom
            final_tin = extract_geom(tin_df)
            
            #final_tin = final_tin.overlay(f) #(final_tin.geometry.centroid.within(f.geometry[0])) |
                        
            triangles = pd.concat([triangles, final_tin], ignore_index=True)

        else:
            print('Using only BFE pts')
            # Concat and Triangulate
            all_pts_multigeom = MultiPoint(bfe_pts_84.geometry.to_list())

            tin = triangulate(all_pts_multigeom)
            tin_df = g(tin, 4326)
        
            # Extract Geom
            final_tin = extract_geom(tin_df)
            #final_tin = final_tin.overlay(f) #(final_tin.geometry.centroid.within(f.geometry[0])) |
                                
            triangles = pd.concat([triangles, final_tin], ignore_index=True)

        
        except:
            print('Problem with Poly: ', i)
            p_er_dict = {i: f['geometry'][0]}
            error = gpd.GeoSeries(p_er_dict, crs=26913)
            poly_errors = pd.concat([poly_errors, error])
            continue

if poly_errors:
    poly_errors = poly_errors.to_crs(4326)
        



print('______________________')
print('Writing Outputs in WGS84: ')
triangles.to_file(f'{out}_Triangles.shp')
print(f'{out}_Triangles')
poly_errors.to_file(f'{out}_polyerrors.shp')
print(f'{out}_polyerrors')
fsp_s.to_file(f'{out}_FSP_split.shp')
print(f'{out}_FSP_split')
bfe.to_file(f'{out}_BFE.shp')
print(f'{out}_BFE')
bfe_ext.to_file(f'{out}_BFE_Ext.shp')
print(f'{out}_BFE_Ext')
print('_______________________')
eta = (time.time() - s) / 60
print('Process Complete:', eta, 'minutes')
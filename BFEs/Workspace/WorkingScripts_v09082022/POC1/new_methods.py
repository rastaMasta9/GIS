#!/usr/bin/env python
# coding: utf-8

import os
import geopandas as gpd
from osgeo import gdal
import numpy as np
import pandas as pd
import math
import itertools
import warnings
warnings.filterwarnings('ignore')
from shapely.geometry import Polygon
from shapely.geometry import MultiPoint
from shapely.geometry import MultiPolygon
from shapely.geometry import MultiLineString
from shapely.geometry import LineString
from shapely.geometry import Point
from shapely.ops import triangulate
from shapely.ops import polygonize
from shapely.ops import transform
from shapely.ops import linemerge
from shapely.ops import unary_union
from shapely.ops import nearest_points
from shapely.ops import split
from centerline.geometry import Centerline
from osmnx.utils_geo import _quadrat_cut_geometry
from shapely.affinity import scale

# 0: Geodataframe creation
def g(geom, crs):
    try:
        df = gpd.GeoDataFrame(geometry=geom, crs=crs)
        return df
    except:
        df = gpd.GeoDataFrame(geometry=[geom], crs=crs)
        return df


# 1: Extend BFES Across FSP POLYGON
# 1-1: 
def getExtrapoledLine(p1,p2):
    """Creates a line extrapoled in p1->p2 direction"""
    EXTRAPOL_RATIO = 10
    a = p1
    b = (p1[0]+EXTRAPOL_RATIO*(p2[0]-p1[0]), p1[1]+EXTRAPOL_RATIO*(p2[1]-p1[1]) )
    return LineString([a,b])

# 1-2: 
def get_ext_linestring(fr_pt, ba_pt):
    """Handles Linestrings that are already extending past fsp poly.
    Produce new line linestring"""
    all_coords = []
    pt_group = [fr_pt, ba_pt]

    for l in pt_group:
        if l.geom_type == 'MultiPoint':
            explode = [x for x in l]
            e_coords = [explode[0].coords[0], explode[1].coords[0]]
        
            all_coords.append(e_coords[0])
            all_coords.append(e_coords[1])
        else:
            coords = l.coords[0]
    
            all_coords.append(list(coords))
    nline = LineString(all_coords)
    return nline

#1-3: Applies 1-1, 1-2 to BFE df
def extend_bfe(bfe, poly):
    bfe_e = bfe.copy()
    poly_exterior = poly.geometry.boundary[0]
    bfe_e['points'] = bfe_e['geometry'].apply(lambda x: list((x.boundary[0].coords[0], x.boundary[1].coords[0])))
    bfe_e['fr_ext'] = bfe_e.apply(lambda x: getExtrapoledLine(x['points'][0], x['points'][1]), axis=1)
    bfe_e['ba_ext'] = bfe_e.apply(lambda x: getExtrapoledLine(x['points'][1], x['points'][0]), axis=1)
    bfe_e['fr_pt'] = bfe_e.apply(lambda x: poly_exterior.intersection(x['fr_ext']), axis=1)
    bfe_e['ba_pt'] = bfe_e.apply(lambda x: poly_exterior.intersection(x['ba_ext']), axis=1)

    
    bfe_e['ext_geom'] = bfe_e.apply(lambda x: get_ext_linestring(x['fr_pt'], x['ba_pt']), axis=1)
    bfe_e['ext_geom_scale'] = bfe_e.apply(lambda x: scale(x.geometry, xfact=1.2, yfact=1.2), axis=1)
    
    bfe_e = bfe_e[['ELEV', 'ext_geom_scale']]
    bfe_e.rename(columns={'ext_geom_scale': 'geometry'}, inplace=True)
    bfe_e.set_geometry('geometry', crs=26913, inplace=True)
    
    return bfe_e
    

# 2: SPLIT FSP BY BFE
def split_fsp(fsp, bfe):
    fsp_line = fsp.geometry.boundary[0]
    fsp_line = g(fsp_line, 26913)
    fsp_line = fsp_line.explode()
    lines = bfe['geometry'].to_list() + fsp_line['geometry'].to_list()
    
    merge = linemerge(lines)
    union = unary_union(merge)
    poly = polygonize(union)
    df = g([p for p in poly], 26913)
    df.reset_index(inplace=True)
    
    return df

# 3 ELEV Field to Z-Geometry
def ELEV_2geom(df):
    """Takes Z geometry from interpolated dataframe and 
    creates attribute table"""

    # ADD Z VALUE FROM ELEV FIELD
    df['geometry'] = df['geometry'].apply(lambda p: transform(lambda x, y: (x, y, df.loc[df['geometry'] == p, 'ELEV']), p))

    return df[['ELEV', 'geometry']]

# 4: Getting Z geom for BFE points
def bfe_zpts(bfe):
    bfe['points'] = bfe.apply(lambda x: [y for y in x['geometry'].coords], axis=1)
    geom = [Point(x,y) for coords in bfe['points'].to_list() for x, y in coords]
    
    bfe_pts = g(geom, 26913)
    bfe_pts_sj = (bfe_pts.sjoin(bfe).pipe(ELEV_2geom))

    return bfe_pts_sj

    
# 5: FSP Simplify
def fsp_pts_simplify(fsp_split, tolerance):
    fsp_mask_sim = list(fsp_split.geometry.simplify(tolerance=tolerance))
    geom = [Point(x,y) for coords in fsp_mask_sim for x, y in coords.exterior.coords]
    fsp_pts = g(geom, 26913)
    
    return fsp_pts

# 6: Point Z interpolation
def IDW(pts, bfe, power, bfe_count):
    # Create centroid field
    bfe['centroid'] = bfe.apply(lambda x: Point([x['geometry'].centroid.coords[0]]), axis=1)

    if bfe_count == 1:
        bfe1_elev = bfe.iloc[0]['ELEV']
        for i, geo in pts.iterrows():
            d1 = bfe.iloc[0].centroid.distance(geo['geometry'])
            
            elev = ((bfe1_elev/d1**power) / (1/d1**power))
            pts.loc[i, 'ELEV'] = elev

    elif bfe_count == 2:
        bfe1_elev = bfe.iloc[0]['ELEV']
        bfe2_elev = bfe.iloc[1]['ELEV']

        for i, geo in pts.iterrows():
            d1 = bfe.iloc[0].centroid.distance(geo['geometry'])
            d2 = bfe.iloc[1].centroid.distance(geo['geometry'])
        
            elev = ((bfe1_elev/d1**power) + (bfe2_elev/d2**power)) \
            / ((1/d1**power) + (1/d2**power))
            
            pts.loc[i, 'ELEV'] = elev

    elif bfe_count == 3:
        bfe1_elev = bfe.iloc[0]['ELEV']
        bfe2_elev = bfe.iloc[1]['ELEV']
        bfe3_elev = bfe.iloc[2]['ELEV']

        for i, geo in pts.iterrows():
            d1 = bfe.iloc[0].centroid.distance(geo['geometry'])
            d2 = bfe.iloc[1].centroid.distance(geo['geometry'])
        
            elev = ((bfe1_elev/d1**power) + (bfe2_elev/d2**power)) \
            / ((1/d1**power) + (1/d2**power))
            
            pts.loc[i, 'ELEV'] = elev
    
    else:
        print('More than 3 BFEs!!')
        
    return pts
        

# 7: Extract geometry into attribute table
def extract_geom(df):
# EXTRACT X, Y, Z PER POINT FOR EACH POLYGON IN RESPECTIVE FIELDS
    for i, r in df.iterrows():
        idx = 0
        for pt in list(r['geometry'].exterior.coords):
            lat = pt[1]
            long = pt[0]
            z = pt[2]
            df.loc[i, f'pt_{idx}_LAT'] = lat
            df.loc[i, f'pt_{idx}_LONG'] = long
            df.loc[i, f'pt_{idx}_Z'] = z
            idx+=1
    
    return df

# 8: MultinLine Removal -- Set each BFE to a single LineString
def remove_multiline_BFE(bfe_set):
    if 'MultiLineString' in bfe_set.geom_type.to_list():
        multi_geom = bfe_set.loc[bfe_set.geom_type == 'MultiLineString']
        lines = multi_geom.explode()
        max_line = lines.loc[lines.geometry.length == max(lines.geometry.length)]
        bfe_set = bfe_set.loc[~bfe_set.index.isin(multi_geom.index)]
        bfe_set = pd.concat([bfe_set, max_line], ignore_index=True)
    else:
        pass
        
    return bfe_set
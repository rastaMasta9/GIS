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


# 1: BFE Extension
# 1-1: Extrapolate extension to Find Intersection Across FSP POLYGON
def getExtrapoledLine(p1,p2):
    """Creates a line extrapoled in p1->p2 direction"""
    EXTRAPOL_RATIO = 5
    a = p1
    b = (p1[0]+EXTRAPOL_RATIO*(p2[0]-p1[0]), p1[1]+EXTRAPOL_RATIO*(p2[1]-p1[1]) )
    return LineString([a,b])

# 1-2: Apply Extension to BFE df and buffer new linstring geom to include additional point just outside FSP poly
def bfe_extend(bfe, fsp):
    bfec = bfe.copy()
    poly_exterior = fsp.geometry.boundary[0]
    bfec['cuts'] = bfec.apply(lambda x: math.ceil(x['geometry'].length / 2), axis=1)
    bfec['points'] = bfec.apply(lambda x: [x.geometry.interpolate(i/x['cuts'], normalized=True) for i in range(1, x['cuts'])], axis=1)

    bfec['fr_ext'] = bfec.apply(lambda x: getExtrapoledLine(x['points'][1].coords[0], x['points'][0].coords[0]), axis=1)
    bfec['ba_ext'] = bfec.apply(lambda x: getExtrapoledLine(x['points'][-2].coords[0], x['points'][-1].coords[0]), axis=1)
    bfec['fr_pt'] = bfec.apply(lambda x: poly_exterior.intersection(x['fr_ext']), axis=1)
    bfec['ba_pt'] = bfec.apply(lambda x: poly_exterior.intersection(x['ba_ext']), axis=1)

    bfec = bfec.loc[(bfec['fr_pt'].geom_type != 'LineString') & (bfec['ba_pt'].geom_type != 'LineString')]

    bfec['fr_pt_near'] = bfec.apply(lambda x: [nearest_points(x['fr_pt'], x['geometry'])[0]], axis=1)
    bfec['ba_pt_near'] = bfec.apply(lambda x: [nearest_points(x['ba_pt'], x['geometry'])[0]], axis=1)
    bfec['new_geom'] = bfec.apply(lambda x: LineString(x['fr_pt_near'] + x['points'] + x['ba_pt_near']), axis=1)
    bfec['new_pts'] = bfec.apply(lambda x: [Point(y) for y in x['new_geom'].coords], axis=1)

    bfec_ex = bfec[['ELEV', 'new_pts', 'geometry']] # At this point Line has been snapped to Polygon. But wait, there's more!

    bfec_ex['ends'] = bfec_ex.apply(lambda x: [Point(x.geometry.boundary[0].coords[0]), Point(x.geometry.boundary[1].coords[0])], axis=1)
    bfec_ex['buff'] = bfec_ex.apply(lambda x: (x['ends'][0].buffer(10), x['ends'][1].buffer(10)), axis=1)
    bfec_ex['fr_diff'] = bfec_ex.apply(lambda x: x['buff'][0].difference(fsp.geometry[0]), axis=1)
    bfec_ex['ba_diff'] = bfec_ex.apply(lambda x: x['buff'][1].difference(fsp.geometry[0]), axis=1)
    bfec_ex['frpt'] = bfec_ex.apply(lambda x: [x['fr_diff'].representative_point()], axis=1)
    bfec_ex['bapt'] = bfec_ex.apply(lambda x: [x['ba_diff'].representative_point()], axis=1)

    bfec_ex['new_geom'] = bfec_ex.apply(lambda x: LineString(x['frpt'] + x['new_pts'] + x['bapt']), axis=1)

    bfec_ex = bfec_ex[['ELEV', 'new_geom']]
    bfec_ex.rename(columns={'new_geom': 'geometry'}, inplace=True)
    bfec_ex = bfec_ex.set_geometry('geometry', crs=26913)

    return bfec_ex


    

# 2: SPLIT FSP BY BFE
def split_fsp(fsp, bfe):
    fsp_line = fsp.geometry.boundary[0]
    fsp_line = g(fsp_line, 26913)
    fsp_line = fsp_line.explode()

    # Clean up BFE_extension (i.e. ensure only LineStrings occur)
    bfe = bfe.loc[bfe.geom_type == 'LineString']

    # Collect
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
def bfe_zpts(df):
    set = df.copy()
    set['ends'] = set.apply(lambda x: list(x.geometry.boundary), axis=1)
    #set['center'] = set.apply(lambda x: [x.geometry.centroid], axis=1)
    #set['points'] = set.apply(lambda x: x['ends'] + x['center'], axis=1)

    geom = [coords for x in set['ends'].to_list() for coords in x]
    pp = g(geom, 26913)
    pp = g(pp.buffer(1), 26913)
    sj = pp.sjoin(df)
    sj['center'] = sj.apply(lambda x: x.geometry.centroid, axis=1)
    sj = sj[['ELEV', 'center']]
    pts = sj.set_geometry('center', crs=26913)
    pts = sj.rename(columns={'center': 'geometry'})
    
    zpts = ELEV_2geom(pts)

    return zpts


    
# 5: FSP Interpolated points along Extent
def fsp_pts(fsp_split, spacing):
    cuts = math.ceil(fsp_split.length / spacing)
    if cuts <= 15:
        cuts = 30

    geom = [fsp_split.geometry.boundary[0].interpolate(i/cuts, normalized=True) for i in range(1, cuts)]
    fsp_pts = g(geom, 26913)
    
    return fsp_pts

# 6: Point Z interpolation
def IDW(pts, bfe, power, buffer=1):
    """An Inverse Distance Weighted Algorithm to interpolate elevations for points
    that had been generated along the fsp extent. Loops through each bfe that intersects
    the current fsp instance. Breaks the IDW equation into fields making it easier to calcualte
    row-wise within a geopdataframe. """

    for b in range(len(bfe)):
        # Calcualtes distance for each point from each BFE. If Point resides along edge of BFE( within 1-meter), distance is is given NUll then dropped
        pts[f'd_to_{b}'] = pts.apply(lambda x: bfe.iloc[b].geometry.distance(x['geometry']) if bfe.iloc[b].geometry.distance(x['geometry']) > buffer else np.nan, axis=1)
        pts = pts.dropna(how='any')
        # The Nominator calcualtion of the IDW algorithm
        pts[f'nom_{b}'] = pts.apply(lambda x: (bfe.iloc[b]['ELEV'] / x[f'd_to_{b}']**power), axis=1)

        # The Denominator calcualtion of the IDW algorithm
        pts[f'dom_{b}'] = pts.apply(lambda x: (1 / x[f'd_to_{b}']**power), axis=1)

    # Capture Nominator and Denominator Columns to later total
    nom_col = pts.filter(regex='nom')
    dom_col = pts.filter(regex='dom')

    # Sum Correspnding Nominator and Denominator columns. Based on the IDW algorithm
    pts['nom_totals'] = pts[nom_col.columns.to_list()].sum(axis=1)
    pts['dom_totals'] = pts[dom_col.columns.to_list()].sum(axis=1)

    # Final Calculation in the IDW algorithm. (Nom / Dom)
    pts['ELEV'] = pts.apply(lambda x: (x['nom_totals'] / x['dom_totals']), axis=1)
    pts['ELEV'] = pts.apply(lambda x: round(x['ELEV'], 2), axis=1)

    return pts[['ELEV', 'geometry']]
        

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
# 8-1 Get Longest Linestring and replace geom if MultiLineString:
def get_longLine(multi_geom):
    lines = multi_geom.explode()
    lines.reset_index(inplace=True)
    lines_gb = lines.groupby('level_0').apply(lambda x: max(x.geometry.length))
    gb_dict = lines_gb.to_dict()
    lines['max_geom'] = lines['level_0'].map(gb_dict)
    longs = lines.loc[lines['geometry'].length == lines['max_geom'], ['level_0', 'ELEV', 'geometry']]
    longs.set_index('level_0', inplace=True)
    longs.index.name = None
    
    return longs

# 8-2 Accomodate for instances when only one is MultiLinestring
def remove_multiline_BFE(bfe_set):
    if 'MultiLineString' in bfe_set.geom_type.to_list():
        multi_geom = bfe_set.loc[bfe_set.geom_type == 'MultiLineString']
        single = bfe_set.loc[bfe_set.geom_type != 'MultiLineString']
    
        if len(single) == 0:
            longs = get_longLine(multi_geom)
            return longs
            
        else:
            
            longs = get_longLine(multi_geom)
            longs = pd.concat([longs, single])
            return longs
    else:
        return bfe_set
        
    
#!/usr/bin/env python
# coding: utf-8

import os
import glob
import geopandas as gpd
from osgeo import gdal
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings('ignore')
from shapely.geometry import Polygon
from shapely.geometry import MultiPoint
from shapely.geometry import MultiPolygon
from shapely.geometry import MultiLineString
from shapely.geometry import Point
from shapely.geometry import mapping
from shapely.ops import triangulate
from shapely.ops import polygonize_full
from shapely.ops import polygonize
from shapely.ops import transform
from shapely.ops import linemerge
from shapely.ops import unary_union
from shapely.ops import nearest_points
from shapely.ops import split
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
def extend_bfe(bfe, poly):
    bfe_e = bfe.copy()
    bfe_e['points'] = bfe_e.apply(lambda x: [y for y in x['geometry'].coords], axis=1)
    
    bfe_e['dist_2poly'] = bfe_e.apply(lambda x: (Point([x['points'][0]]).distance(poly.geometry[0]),
                                       Point([x['points'][1]]).distance(poly.geometry[0])), axis=1)
    
    bfe_e['geometry'] = bfe_e.apply(lambda x: scale(x.geometry, xfact=1.5, yfact=1.5)
                            if int(x['dist_2poly'][0]) or int(x['dist_2poly'][1]) < 10
                            else (scale(x.geometry, xfact=3, yfact=3)
                                  if int(x['dist_2poly'][0]) | int(x['dist_2poly'][1]) < 40
                            else scale(x.geometry, xfact=5, yfact=5)), axis=1)
    bfe_e['points'] = bfe_e.apply(lambda x: [y for y in x['geometry'].coords], axis=1)
    
    bfe_e.set_geometry('geometry', inplace=True)
    bfe_e = bfe_e[['ELEV', 'geometry']]
    
    return bfe_e
    

# 2: SPLIT FSP BY BFE
def split_fsp(fsp, bfe):
    fsp_line = fsp.geometry.boundary[0]
    fsp_line = g(fsp_line, 26913)
    lines = bfe['geometry'].to_list() + fsp_line['geometry'].to_list()
    
    merge = linemerge(lines)
    union = unary_union(merge)
    poly = polygonize(union)
    df = g([p for p in poly], 26913)
    df.reset_index(inplace=True)
    
    return df

# 3: Split Flowline to Iterate on
def split_flowline(fl, bfe):
    fl_line = fl.geometry.unary_union
    bfe_line = bfe.geometry.unary_union
    
    fl_split = split(linemerge(fl_line), bfe_line)
    fl_df = g([f for f in fl_split], 26913)
    
    return fl_df

# 4: Buffer line segs by approx 1m and then spatial join with original flowline df
def buff_fl_intersection(fl, buffer):
    buff_geom = fl['geometry'].buffer(buffer)
    fl_buff = g(buff_geom, 26913)
    fl_buff.reset_index(inplace=True)

    fl_buff_sj = fl_buff.sjoin(fl)
    fl_buff_sj = fl_buff_sj.loc[fl_buff_sj['index'] != fl_buff_sj['index_right']]
    
    return fl_buff, fl_buff_sj

# 5 Create Adjacent df to parse Forks and Normal Flowline segments.
# Iterate on this df during build
def create_adjacent_gpd(fl_buff_sj):
    adj = fl_buff_sj.groupby('index')['index_right'].apply(list)
    adj = pd.DataFrame(adj)
    adj.reset_index(inplace=True)

    return adj


# 6: Find Forks from Adj flowline segments
def find_forks(df):
    """USE THREES DATAFRAME. Takes flowlines with more than 2 intersections and organizes intersectional list based on
    which indices occur within them. Returns a list of indices that create each fork"""
    
    # Locate intesections where indices occur and sort each list to create flowline indice lists.
    fdict = {}
    for i in df['index']:
        mask = df['index_right'].apply(lambda x: i in x)
        fls = df.loc[mask, 'index'].to_list()
        fls.append(i)
        fls.sort()
        fdict[i] = fls
    
    # Build df and drop simliar lists (i.e. leaves only unqiue lists)
    fdf = pd.DataFrame.from_dict(fdict, orient='index')
    fdf = fdf.drop_duplicates()
    
    # Retrieve each indice found in the df column and place in correspnding list
    forks = []
    for i, r in fdf.iterrows():
        l = [r[0], r[1], r[2]]
        forks.append(l)
    return forks


# 7: Getting Z geom for BFE points
def bfe_zpts(bfe):
    bfe['points'] = bfe.apply(lambda x: [y for y in x['geometry'].coords], axis=1)
    geom = [Point(x,y) for coords in bfe['points'].to_list() for x, y in coords]
    
    bfe_pts = g(geom, 26913)
    bfe_pts_sj = bfe_pts.sjoin(bfe)
    
    # ADD Z VALUE FROM ELEV FIELD
    bfe_pts_sj['geometry'] = bfe_pts_sj['geometry'].apply(
            lambda p: 
            transform(lambda x, y: 
            (x, y, bfe_pts_sj.loc[bfe_pts_sj['geometry'] == p, 'ELEV']), p))
    
    return bfe_pts_sj

# 7: BFE Centroid and interpolate pts along Flowline
def flowline_interpolation(bfe, fl, divisions):
    # get distance between BFEs
    tot_d = bfe.iloc[0].geometry.distance(bfe.iloc[1].geometry)
    cuts = round(tot_d/divisions)
    
    # get bfe centroid and interpolate points with elev
    bfe['centroid'] = bfe.apply(lambda x: Point([x['geometry'].centroid.coords[0]]), axis=1)
    bfe_center = bfe.set_geometry('centroid')
    bfe_center = bfe_center[['ELEV', 'centroid']]
    
    # interpolate points at fixed distance 'cuts'
    splitter = MultiPoint([fl.geometry.interpolate(i/cuts, normalized=True) for i in range(1, cuts)])
    
    interp_geom = [s for s in splitter]
    interp_df = g(interp_geom, 26913)
    
    forward = (0, 1)
    reverse = (1, 0)
    orient = [forward, reverse]
    for o in orient:
        bfe1_elev = bfe.iloc[o[0]]['ELEV']
        bfe2_elev = bfe.iloc[o[1]]['ELEV']
        dist = pd.Series()
        
        for i in interp_df['geometry']:
            d = bfe.iloc[o[0]]['centroid'].distance(i)
            dist = pd.concat([dist, pd.Series(d)], ignore_index=True)
        total = dist.sum()

        for i, geo in interp_df.iterrows():
            d = bfe.iloc[o[0]]['centroid'].distance(geo['geometry'])
            elev = ((bfe1_elev/d) + (bfe2_elev/(total - d))) / ((1/d) + (1/(total - d)))
            interp_df.loc[i, f'ELEV_{o[0]}'] = elev
        
    interp_df['ELEV'] = interp_df.apply(lambda x: ((x['ELEV_0'] + x['ELEV_1'])/2), axis=1)
    interp_df.reset_index(inplace=True)
    
    # ADD Z VALUE FROM ELEV FIELD
    interp_df['geometry'] = interp_df['geometry'].apply(
            lambda p: 
            transform(lambda x, y: 
            (x, y, interp_df.loc[interp_df['geometry'] == p, 'ELEV']), p))
    
    return interp_df[['index', 'ELEV', 'geometry']]
    
    
# 7: FSP Simplify
def fsp_pts_simplify(fsp_split, fl_interpolate_POINT, tolerance):
    fl_pt = g(fl_interpolate_POINT.iloc[0].geometry, 26913)

    fsp_mask = fsp_split.sjoin(fl_pt)
    fsp_mask_sim = list(fsp_mask.geometry.simplify(tolerance=tolerance))
    geom = [Point(x,y) for coords in fsp_mask_sim for x, y in coords.exterior.coords]
    
    fsp_pts = g(geom, 26913)
    
    return fsp_pts

# 8: FSP Point Z interpolation
def IDW(bfe, pts):
    forward = (0, 1)
    reverse = (1, 0)
    orient = [forward, reverse]
    for o in orient:
        bfe1_elev = bfe.iloc[o[0]]['ELEV']
        bfe2_elev = bfe.iloc[o[1]]['ELEV']
        dist = pd.Series()

        for i in pts['geometry']:
            d = bfe.iloc[o[0]]['centroid'].distance(i)
            dist = pd.concat([dist, pd.Series(d)], ignore_index=True)
        total = dist.sum()

        for i, geo in pts.iterrows():
            d = bfe.iloc[o[0]]['centroid'].distance(geo['geometry'])
            elev = ((bfe1_elev/d) + (bfe2_elev/(total - d))) / ((1/d) + (1/(total - d)))
            pts.loc[i, f'ELEV_{o[0]}'] = elev

    pts['ELEV'] = pts.apply(lambda x: ((x['ELEV_0'] + x['ELEV_1'])/2), axis=1)
    pts['ELEV'] = pts['ELEV'].fillna(0)
    pts.reset_index(inplace=True)

    # ADD Z VALUE FROM ELEV FIELD
    pts['geometry'] = pts['geometry'].apply(
            lambda p: 
            transform(lambda x, y: 
            (x, y, pts.loc[pts['geometry'] == p, 'ELEV']), p))

    return pts[['index', 'ELEV', 'geometry']]

# 9: Extract geometry into attribute table
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
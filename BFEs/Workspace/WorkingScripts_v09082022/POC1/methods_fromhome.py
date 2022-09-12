#!/usr/bin/env python
# coding: utf-8

import os
import geopandas as gpd
from osgeo import gdal
import numpy as np
import pandas as pd
import itertools
import warnings
warnings.filterwarnings('ignore')
from shapely.geometry import Polygon
from shapely.geometry import MultiPoint
from shapely.geometry import MultiPolygon
from shapely.geometry import MultiLineString
from shapely.geometry import Point
from shapely.ops import triangulate
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
                            if int(x['dist_2poly'][0]) | int(x['dist_2poly'][1]) < 10
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
    bfe_l = bfe.copy()
    fl_line = fl.geometry.unary_union
    bfe_line = bfe_l.geometry.unary_union
    
    fl_split = split(linemerge(fl_line), bfe_line)
    fl_df = g([f for f in fl_split], 26913)
    
    return fl_df

# 4: Find Forks from Adj flowline segments
def find_forks(flowline):
    forks_list = []
    for i, f in flowline.iterrows():
        p1, p2 = f.geometry.boundary
        p1b = g(p1.buffer(5), 26913)
        p2b = g(p2.buffer(5), 26913)
        p1j = p1b.sjoin(flowline)
        p2j = p2b.sjoin(flowline)

        if p1j.shape[0] == 3:
            forks_list.append(p1j['index_right'].to_list())
        if p2j.shape[0] == 3:
            forks_list.append(p2j['index_right'].to_list())
    
    forks_list.sort()
    forks_list = list(f for f,_ in itertools.groupby(forks_list))
    
    return forks_list

# 5: Getting Z geom for BFE points
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

# 6: Create points from non-fork LineStrings
def interp_pts_fromLine(bfe, line_union, divisions):
# get distance between BFEs
    tot_d = bfe.iloc[0].geometry.distance(bfe.iloc[1].geometry)
    cuts = round(tot_d/divisions)

    # interpolate points at fixed distance 'cuts'
    splitter = MultiPoint([line_union.geometry.interpolate(i/cuts, normalized=True) for i in range(1, cuts)])

    interp_geom = [s for s in splitter]
    interp_pts = g(interp_geom, 26913)

    return interp_pts

# 7: Create points from Fork or Fake Fork Linestrings
def interp_pts_fromFork(bfe, 
                        line_union, 
                        line_segs, 
                        buff_segments, 
                        divisions, 
                        non_fork):
    # get distance between BFEs
    if non_fork:
        bsegs_sj = buff_segments.sjoin(bfe, how='left')
        ids = list(bsegs_sj[~bsegs_sj.isna().any(axis=1)].index)
        line_segs = line_segs.loc[line_segs.index.intersection(ids)]
        nline_union = g(MultiLineString(line_segs.geometry.to_list()), 26913)

        tot_d = bfe.iloc[0].geometry.distance(bfe.iloc[1].geometry)
        cuts = round(tot_d/divisions)
        # interpolate points at fixed distance 'cuts'
        splitter = [list(nline_union.geometry.interpolate(i/cuts, normalized=True)) for i in range(1, cuts)]

        interp_geom = MultiPoint([s for x in splitter for s in x])
        interp_pts = g(interp_geom, 26913)
        interp_pts = interp_pts.explode()
    
        return interp_pts

    else:
        tot_d1 = bfe.iloc[0].geometry.distance(bfe.iloc[1].geometry)
        tot_d2 = bfe.iloc[0].geometry.distance(bfe.iloc[2].geometry)
        tot_d3 = bfe.iloc[1].geometry.distance(bfe.iloc[2].geometry)
        d_dict = {'1': tot_d1, '2': tot_d2, '3': tot_d3}

        max_fork_dist = max(d_dict.values())
        cuts = round(max_fork_dist/divisions)

        # interpolate points at fixed distance 'cuts'
        splitter = [list(line_union.geometry.interpolate(i/cuts, normalized=True)) for i in range(1, cuts)]

        interp_geom = MultiPoint([s for x in splitter for s in x])
        interp_pts = g(interp_geom, 26913)
        interp_pts = interp_pts.explode()

        return interp_pts

# 8: BFE Centroid and interpolate pts along Flowline
def flowline_interpolation(bfe, 
                            line_union, 
                            line_segs, 
                            buff_segments, 
                            divisions, 
                            power, 
                            fork='no', 
                            non_fork=False):
    if fork == 'yes':
        interp_f_pts = interp_pts_fromFork(bfe, line_union, line_segs, buff_segments, divisions, non_fork)
        interp_f_df = IDW_Forks(bfe, interp_f_pts, power, non_fork)
        return interp_f_df

    else:
        interp_pts = interp_pts_fromLine(bfe, line_union, divisions)
        interp_df = IDW(bfe, interp_pts, power)
        return interp_df

    
# 9: FSP Simplify
def fsp_pts_simplify(fsp_split, fl_interpolate_POINT, tolerance):
    for i in fl_interpolate_POINT['geometry']:
        fl_pt = g(i, 26913)

        fsp_mask = fsp_split.sjoin(fl_pt)
        if fsp_mask.shape[0] != 0:
            fsp_mask_sim = list(fsp_mask.geometry.simplify(tolerance=tolerance))
            geom = [Point(x,y) for coords in fsp_mask_sim for x, y in coords.exterior.coords]
            
            fsp_pts = g(geom, 26913)
            break
        
    return fsp_pts

# 10: Normal 2 BFE Point Z interpolation
def IDW(bfe, pts, power):
    # Create centroid field
    bfe['centroid'] = bfe.apply(lambda x: Point([x['geometry'].centroid.coords[0]]), axis=1)

    bfe1_elev = bfe.iloc[0]['ELEV']
    bfe2_elev = bfe.iloc[1]['ELEV']
    
    for i, geo in pts.iterrows():
        d1 = bfe.iloc[0].centroid.distance(geo['geometry'])
        d2 = bfe.iloc[1].centroid.distance(geo['geometry'])
        
        elev = ((bfe1_elev/d1**power) + (bfe2_elev/d2**power)) \
        / ((1/d1**power) + (1/d2**power))
        pts.loc[i, 'ELEV'] = elev

    # ADD Z VALUE FROM ELEV FIELD
    pts['geometry'] = pts['geometry'].apply(
            lambda p: 
            transform(lambda x, y: 
            (x, y, pts.loc[pts['geometry'] == p, 'ELEV']), p))

    return pts[['ELEV', 'geometry']]

#  11: Fork BFE Point Z Interpolation
def IDW_Forks(bfe, interp_pts, power, non_fork):
    if non_fork:
        interp_pts = IDW(bfe, interp_pts, power)
    
    else:
        # Create centroid field
        bfe['centroid'] = bfe.apply(lambda x: Point([x['geometry'].centroid.coords[0]]), axis=1)

        bfe1_elev = bfe.iloc[0]['ELEV']
        bfe2_elev = bfe.iloc[1]['ELEV']
        bfe3_elev = bfe.iloc[2]['ELEV']

        for i, geo in interp_pts.iterrows():
            d1 = bfe.iloc[0].centroid.distance(geo['geometry'])
            d2 = bfe.iloc[1].centroid.distance(geo['geometry'])
            d3 = bfe.iloc[2].centroid.distance(geo['geometry'])
            
            elev = ((bfe1_elev/d1**power) + (bfe2_elev/d2**power) + (bfe3_elev/d3**power)) \
            / ((1/d1**power) + (1/d2**power) + (1/d3**power))
            interp_pts.loc[i, 'ELEV'] = elev

        # ADD Z VALUE FROM ELEV FIELD
        interp_pts['geometry'] = interp_pts['geometry'].apply(
                lambda p: 
                transform(lambda x, y: 
                (x, y, interp_pts.loc[interp_pts['geometry'] == p, 'ELEV']), p))

    return interp_pts[['ELEV', 'geometry']]

# 12: Extract geometry into attribute table
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

# 13: Union Forks (3 segs. to 1 = Union)
def union_fork(ids, fl, fl_buff):
    fork_segs = fl.loc[fl.index.intersection(ids)]
    fork_segs_buff = fl_buff.loc[fl_buff.index.intersection(ids)]
    fork_union = g(MultiLineString(fork_segs.geometry.to_list()), 26913)
    fork_buff_union = g([fork_segs_buff.geometry.unary_union], 26913)
    
    return fork_segs, fork_union, fork_segs_buff, fork_buff_union

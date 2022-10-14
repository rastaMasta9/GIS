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

    # Create points along each LineString to gain x, y for LineString shape
    bfec['cuts'] = bfec.apply(lambda x: math.ceil(x['geometry'].length / 2), axis=1)
    bfec['points'] = bfec.apply(lambda x: [x.geometry.interpolate(i/x['cuts'], normalized=True) for i in range(1, x['cuts'])], axis=1)

    # Tiny BFEs should just be removed (< 1-meter)
    bfec['smalls'] = bfec.apply(lambda x: 1 if len(x['points']) < 3 else np.nan, axis=1)
    bfec = bfec.loc[bfec['smalls'].isnull()]
    
    # Initiate extrapolation on each 'end' of the LineString, extending in the outward direction
    bfec['fr_ext'] = bfec.apply(lambda x: getExtrapoledLine(x['points'][1].coords[0], x['points'][0].coords[0]), axis=1)
    bfec['ba_ext'] = bfec.apply(lambda x: getExtrapoledLine(x['points'][-2].coords[0], x['points'][-1].coords[0]), axis=1)

    # Get Point geometry where the Extrapolated Line eventually intersects the FSP polygon
    bfec['fr_pt'] = bfec.apply(lambda x: poly_exterior.intersection(x['fr_ext']), axis=1)
    bfec['ba_pt'] = bfec.apply(lambda x: poly_exterior.intersection(x['ba_ext']), axis=1)

    # Ensure that the point is either a POINT or MULTIPOINT.
    """ It may at times return a MULTIPOINT because the extrapolated line could interesect
    multiple parts of the FSP Polygon"""
    bfec = bfec.loc[(bfec['fr_pt'].geom_type != 'LineString') & (bfec['ba_pt'].geom_type != 'LineString')]

    # Get POINT that is closest to the original line string (i.e. this accomdates for the MULTIPOINT geometry)
    bfec['fr_pt_near'] = bfec.apply(lambda x: [nearest_points(x['fr_pt'], x['geometry'])[0]], axis=1)
    bfec['ba_pt_near'] = bfec.apply(lambda x: [nearest_points(x['ba_pt'], x['geometry'])[0]], axis=1)

    # Create new Line String from intersection points and the original line string geometry
    bfec['new_geom'] = bfec.apply(lambda x: LineString(x['fr_pt_near'] + x['points'] + x['ba_pt_near']), axis=1)
    bfec['new_pts'] = bfec.apply(lambda x: [Point(y) for y in x['new_geom'].coords], axis=1)
    bfec = bfec[['ELEV', 'new_pts', 'new_geom']]
    bfec.rename(columns={'new_geom': 'geometry'}, inplace=True)
    
    bfec_ex = bfec.set_geometry('geometry', crs=26913) # At this point Line has been snapped to Polygon. But wait, there's more!

    # Get End point geometry
    bfec_ex['ends'] = bfec_ex.apply(lambda x: x.geometry.boundary, axis=1)
    bfec_ex = bfec_ex.loc[~bfec_ex['ends'].is_empty] # This is an anomaly. Needs Investigating. Some geometries are empty. WTF

    # Pull MULTIPOINT geom into seprated Points as tuple
    bfec_ex['ends_pts'] = bfec_ex.apply(lambda x: [Point(x['ends'][0].coords[0]), Point(x['ends'][1].coords[0])], axis=1)

    # Buffer each Point
    bfec_ex['buff'] = bfec_ex.apply(lambda x: (x['ends_pts'][0].buffer(1), x['ends_pts'][1].buffer(1)), axis=1)

    # Return the geometric difference between the intersection of the buffer and FSP_Polygon
    bfec_ex['fr_diff'] = bfec_ex.apply(lambda x: x['buff'][0].difference(fsp.geometry[0]), axis=1)
    bfec_ex['ba_diff'] = bfec_ex.apply(lambda x: x['buff'][1].difference(fsp.geometry[0]), axis=1)

    # Create a random point that is garunteed to lie within the new buffered zone
    bfec_ex['frpt'] = bfec_ex.apply(lambda x: [x['fr_diff'].representative_point()], axis=1)
    bfec_ex['bapt'] = bfec_ex.apply(lambda x: [x['ba_diff'].representative_point()], axis=1)

    # Add rep. point to the new geometry. This adds the slight extenion over the FSP Polygon. Hoorah!
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
    """Takes value in ELEV field and add it to geometry """
    # Reproject to WGS84

    # ADD Z VALUE FROM ELEV FIELD
    df['geometry'] = df['geometry'].apply(lambda p: transform(lambda x, y: (x, y, df.loc[df['geometry'] == p, 'ELEV']), p))

    return df[['ELEV', 'geometry']]

# 4: Getting Z geom for BFE points
def bfe_zpts(df):
    set = df.copy()
    set['ends'] = set.apply(lambda x: list(x.geometry.boundary), axis=1)
    
    # Get centroid
    set['center'] = set.apply(lambda x: x.geometry.representative_point(), axis=1)

    # Create df from end points of each BFE
    ends = [coords for x in set['ends'].to_list() for coords in x]
    geom = ends + set['center'].to_list()
    ppp = g(geom, 26913)

    # This is stupid
    # Required to ensure you will always hit the correspnding BFE during the spatial join
    pp = g(ppp.buffer(1), 26913)
    sj = pp.sjoin(df)
    sj['center'] = sj.apply(lambda x: x.geometry.centroid, axis=1)

    sj = sj[['ELEV', 'center']]
    
    pts = sj.set_geometry('center', crs=26913)
    ptsutm = pts.rename(columns={'center': 'geometry'})
    ptsutm = ptsutm[['ELEV', 'geometry']]

    # Reproject back to WGS 84
    pts84 = pts.to_crs(4326)
    pts84 = pts84.rename(columns={'center': 'geometry'})
    pts84 = pts84[['ELEV', 'geometry']]
    
    zpts84 = ELEV_2geom(pts84)
    zptsutm = ELEV_2geom(ptsutm) # See Function Above

    return zptsutm, zpts84


    
# 5: FSP Interpolated points along Extent
def fsp_pts(fsp_split, bfe_pts, bfe_set, diff_area):
    fsp_pts = None
    bfe_tri = g(triangulate(MultiPoint(bfe_pts.geometry.to_list())), 26913)
    bfe_poly = g(bfe_tri.unary_union, 26913)
    f_diff = fsp_split.overlay(bfe_poly, how='difference')
    diff = f_diff.explode()
    
    if diff.area.max() > diff_area:
        points = 8
        f_line_diff = g(fsp_split.geometry.boundary, 26913).overlay(bfe_set, how='difference')
        geom = [f_line_diff.geometry[0].interpolate(i/points, normalized=True) for i in range(1, points)]
        fsp_pts = g(geom, 26913)

    return fsp_pts

"""
    
# 5: FSP Interpolated points along Extent
def fsp_pts(fsp_split, spacing):
 Spacing parameter is currently set to 50 (i.e. a point every 50 meters)
    Accomodation for smaller Polygons resets the cut variable to a higher Cut rate 
    to ensure more points along smaller extent. Otherwise, they may get removed during the buffer
    of the IDW algorithm
    cuts = math.ceil(fsp_split.length / spacing)
    if 1 < cuts < 5: # poly length between 150 and 500 meters. Double Cuts 
        cuts = cuts*2
    if cuts <= 1:
        cuts = cuts*8 # poly length less than 150, quadrupile.

    geom = [fsp_split.geometry.boundary[0].interpolate(i/cuts, normalized=True) for i in range(1, cuts)]
    fsp_pts = g(geom, 26913)
    
    return fsp_pts
"""

# 6: Point Z interpolation
def IDW(pts, bfe, power, buffer=1):
    """ An Inverse Distance Weighted Algorithm to interpolate elevations for points
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
    pts84 = pts.to_crs(4326)

    return pts84[['ELEV', 'geometry']]
        

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
    # Flatten MultiLine into indivdual LineStrings
    lines = multi_geom.explode()
    lines.reset_index(inplace=True)

    # Groupby BFE index and its Longest indivdual LineString
    lines_gb = lines.groupby('level_0').apply(lambda x: max(x.geometry.length))
    gb_dict = lines_gb.to_dict()

    # Map Value back to original df
    lines['max_geom'] = lines['level_0'].map(gb_dict)

    # Filter Out the other LineStrings. Leave only the Max for each BFE
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
        
    
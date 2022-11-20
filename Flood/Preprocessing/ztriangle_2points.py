#!/usr/bin/env python
# coding: utf-8

# In[2]:


import geopandas as gpd
import pandas as pd

df = gpd.read_file(input('Enter Z-Triangle Shapefile: '))
out = input('Enter Output Shapefile: ')
crs = df.crs.to_epsg()

df0 = df[['pt_0_LAT', 'pt_0_LONG', 'pt_0_Z']]
df1 = df[['pt_1_LAT', 'pt_1_LONG', 'pt_1_Z']]
df2 = df[['pt_2_LAT', 'pt_2_LONG', 'pt_2_Z']]

pt0 = gpd.GeoDataFrame(df0, geometry=gpd.points_from_xy(df0['pt_0_LONG'], df0['pt_0_LAT']), crs=crs)[['pt_0_Z', 'geometry']]
pt1 = gpd.GeoDataFrame(df1, geometry=gpd.points_from_xy(df1['pt_1_LONG'], df1['pt_1_LAT']), crs=crs)[['pt_1_Z', 'geometry']]
pt2 = gpd.GeoDataFrame(df2, geometry=gpd.points_from_xy(df2['pt_2_LONG'], df2['pt_2_LAT']), crs=crs)[['pt_2_Z', 'geometry']]

pt0.rename(columns={'pt_0_Z': 'ELEV'}, inplace=True)
pt1.rename(columns={'pt_1_Z': 'ELEV'}, inplace=True)
pt2.rename(columns={'pt_2_Z': 'ELEV'}, inplace=True)

pts = pd.concat([pt0, pt1, pt2], ignore_index=True)
pts.to_file(out)


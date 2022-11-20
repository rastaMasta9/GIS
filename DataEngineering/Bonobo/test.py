import bonobo
import json
import requests
import pandas as pd
import geopandas as gpd
#from test_extract import *



def extract():
    # Extract from National Park API
    apikey = 'GMs4xrz7WjDCKhiknJtfa6APYdwHawcca4sf6uWE'

    # End point for campgrounds in RockyMountain National Park
    api_url = 'https://developer.nps.gov/api/v1/campgrounds?'

    req = requests.get(api_url, params={
        'api_key': apikey, 
        'stateCode': '08',
        'parkCode': 'romo',
        'limit': 10
        }
    )
    res = req.json().get('data')

    return res

def transform(args):
    # Transform Data into geodataframe
    cg = pd.DataFrame()
    ls = []
    for r in args:
        data = [r['name'], r['latitude'], r['longitude']]
        ls.append(data)

    df = pd.DataFrame(ls, columns=['Campground', 'Lat', 'Long'])
    cg = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df['Long'], df['Lat']), crs=4326)
    
    # Reproject to Web Mercator for use in web map
    cg = cg.to_crs(3857)
    cg = cg[['Campground', 'geometry']]

    yield cg

# Load/Write to File
def load(cg):
    yield cg.to_file('Rocky_Mountain_Campgrounds.shp')

def get_graph():
    graph = bonobo.Graph()
    graph.add_chain(
        extract,
        transform,
        load,
        bonobo.Limit(10)
    )
    return graph

def get_services(**options):
    """
    This function builds the services dictionary, which is a simple dict of names-to-implementation used by bonobo
    for runtime injection.

    It will be used on top of the defaults provided by bonobo (fs, http, ...). You can override those defaults, or just
    let the framework define them. You can also define your own services and naming is up to you.

    :return: dict
    """
    return {}

if __name__ == '__main__':
    parser = bonobo.get_argument_parser()
    with bonobo.parse_args(parser) as options:
        bonobo.run(get_graph(**options),
        services=get_services(**options)
    )


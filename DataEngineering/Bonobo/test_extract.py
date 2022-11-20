# Create transformations
def extract():
    import requests
    import json
    

    # datastored in a list object where each campground is a dictionary
    res = req.json().get('data')

    return res
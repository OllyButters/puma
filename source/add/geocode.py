#! /usr/bin/env python

import csv
import re
import logging
# import json


# Have a go at geocoding the cleaned institute names
# I would expect there to be a lat long for all of them
def geocode(papers):

    print 'Geocoding'

    # Read in lat_long lookup file. Should have an inst name and lat long in each row
    geocode = {}
    lat_long = []
    with open('../config/lat_longs.csv', 'rb') as csvfile:
        f = csv.reader(csvfile)
        for row in f:
            try:
                # Check it is not a comment string first.
                if(re.match('#', row[0])):
                    continue

                # If there is a second element in this row then carry on
                lat_long = {'lat': row[1], 'long': row[2]}
                geocode[row[0]] = lat_long
            except:
                logging.warn('Something went wrong with looking up the lat/long for '+row[0])

    # Actually do the geocoding
    for this_paper in papers:

        try:
            this_paper['Extras']['LatLong'] = geocode[this_paper['Extras']['CleanInstitute']]

        except:
            try:
                print 'Did not find a lat-long for '+this_paper['IDs']['hash']+' '+this_paper['Extras']['CleanInstitute']
                logging.warn('Did not find a lat-long for '+this_paper['IDs']['hash']+' '+this_paper['Extras']['CleanInstitute'])
            except:
                print 'Did not find a lat-long for '+this_paper['IDs']['hash']
                logging.warn('Did not find a lat-long for '+this_paper['IDs']['hash'])

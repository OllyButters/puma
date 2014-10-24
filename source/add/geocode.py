#! /usr/bin/env python


#Have a go at geocoding the cleaned institute names
def geocode(pmids, papers):
    import csv
    import re
    import logging

    #Read in config file
    geocode = {}
    lat_long = []
    with open('../config/lat_longs.csv', 'rb') as csvfile:
        f = csv.reader(csvfile)
        for row in f:
            try:
                #Check it is not a comment string first.
                if(re.match('#',row[0])):
                    continue

                #If there is a second element in this row then carry on
                lat_long = {'lat':row[1], 'long':row[2]}
                geocode[row[0]] = lat_long
            except:
                pass
    #print geocode

    for this_pmid in pmids:
        #print geocode[papers[this_pmid]['Extras']['CleanInstitute']]
        try:
            papers[this_pmid]['Extras']['LatLong'] = geocode[papers[this_pmid]['Extras']['CleanInstitute']]
        except:
            pass

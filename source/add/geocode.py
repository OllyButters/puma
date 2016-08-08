#! /usr/bin/env python

# import csv
# import re
import logging
import os.path
import json
import urllib2
import requests

import config.config as config


# Have a go at geocoding the cleaned institute names
# I would expect there to be a lat long for all of them
def geocode(papers, error_log):

    print 'Geocoding'

    locations_found = 0
    number_done = 1
    for this_paper in papers:

        try:

            # Check Cache
            clean = this_paper['Extras']['CleanInstitute']
            if os.path.isfile(config.cache_dir + "/geodata/" + clean):
                # Load cached data

                cache_file = open(config.cache_dir + "/geodata/" + clean, "r")
                data = cache_file.read()
                split = data.split("#")

                this_paper['Extras']['LatLong'] = {'lat': split[0], 'long': split[1]}
                this_paper['Extras']['country_code'] = split[2]
                this_paper['Extras']['postal_town'] = split[3]

                locations_found += 1

            else:

                # === Look up clean institute geo location ===
                # Form query and get data
                query = 'PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#> SELECT ?item WHERE { ?item rdfs:label "' + this_paper['Extras']['CleanInstitute'] + '"@en }'
                url = 'https://query.wikidata.org/bigdata/namespace/wdq/sparql'

                try:
                    data = requests.get(url, params={'query': query, 'format': 'json'}).json()

                    try:

                        # Process returned JSON to get entity id
                        item_uri = data['results']['bindings'][0]['item']['value']
                        item_uri_components = item_uri.split("/")
                        item_id = item_uri_components[len(item_uri_components)-1]
                        # print( "WD Item id: " + item_id )
                        this_paper['wikidata_item_id'] = item_id

                        # -- USE WIKIDATA ID TO GET GEO-COORDS
                        found_coords = False
                        try:

                            retur = json.load(urllib2.urlopen('https://www.wikidata.org/w/api.php?action=wbgetentities&ids=' + item_id + '&format=json'))
                            p_lon = retur['entities'][item_id]['claims']['P625'][0]['mainsnak']['datavalue']['value']['longitude']
                            p_lat = retur['entities'][item_id]['claims']['P625'][0]['mainsnak']['datavalue']['value']['latitude']

                            locations_found += 1
                            this_paper['Extras']['LatLong'] = {'lat': str(p_lat), 'long': str(p_lon)}

                            found_coords = True

                        except:
                            try:

                                # Check headquaters location
                                p_lon = retur['entities'][item_id]['claims']['P159'][0]['qualifiers']['P625'][0]['datavalue']['value']['longitude']
                                p_lat = retur['entities'][item_id]['claims']['P159'][0]['qualifiers']['P625'][0]['datavalue']['value']['latitude']

                                # print "Location Found from HQ " + papers[this_pmid]['Extras']['CleanInstitute']
                                locations_found += 1
                                this_paper['Extras']['LatLong'] = {'lat': str(p_lat), 'long': str(p_lon)}

                                found_coords = True

                            except:
                                print 'Unable to get geo-data (No HQ P625 Or P625) ' + this_paper['Extras']['CleanInstitute'] + " " + item_id + " (" + str(number_done) + "/" + str(len(papers)) + ")"

                        if found_coords:
                            # Get country for heatmap
                            retur = json.load(urllib2.urlopen('https://maps.googleapis.com/maps/api/geocode/json?latlng=' + this_paper['Extras']['LatLong']['lat']  + ',' + this_paper['Extras']['LatLong']['long']  + '&key=AIzaSyA63o6tsqqAhAB_iPR7foPHEmAU5HMiLe4'))

                            try:
                                comps = retur['results'][0]['address_components']
                                country_short = ""
                                postal_town = ""
                                for comp in comps:
                                    if comp['types'][0] == "country":
                                        country_short = comp['long_name']
                                        this_paper['Extras']['country_code'] = country_short

                                    if comp['types'][0] == "postal_town":
                                        postal_town = comp['long_name']
                                        this_paper['Extras']['postal_town'] = postal_town

                                # Cache Data
                                cache_file = open(config.cache_dir + "/geodata/" + clean, "w")
                                cache_file.write(this_paper['Extras']['LatLong']['lat'] + "#" + this_paper['Extras']['LatLong']['long'] + "#" + country_short + "#" + postal_town)
                                cache_file.close()
                            except:
                                print 'Unable to get geo-data (Google API Quota Reached) ' + this_paper['Extras']['CleanInstitute'] + " (" + str(number_done) + "/" + str(len(papers)) + ")"
                                error_log.logError( this_paper['Extras']['CleanInstitute'] + " Google API Quota Reached!")
                    except:
                        print 'Unable to get geo-data (Probably not on Wikidata) ' + this_paper['Extras']['CleanInstitute'] + " (" + str(number_done) + "/" + str(len(papers)) + ")"
                        error_log.logWarning("Insititue " + this_paper['Extras']['CleanInstitute'] + " not on Wikidata")
                except:
                    print 'Unable to get geo-data (Wikidata Query Failed) ' + this_paper['Extras']['CleanInstitute'] + " (" + str(number_done) + "/" + str(len(papers)) + ")"
                    error_log.logWarning("Wikidata query failed for " + this_paper['Extras']['CleanInstitute'] )

                # === End Look up ===

        except:
            print 'No Clean Institute for ' + this_paper['IDs']['hash'] + " (" + str(number_done) + "/" + str(len(papers)) + ")"
            error_log.logError("Clean Institute Missing for " +  this_paper['IDs']['hash'] )

        number_done += 1

    print "locations found: " + str(locations_found) + "/" + str(number_done)


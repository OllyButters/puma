#! /usr/bin/env python

import csv
# import re
import os.path
import json
import urllib2
import requests

import config.config as config


# Have a go at geocoding the cleaned institute names
# I would expect there to be a lat long for all of them
def geocode(papers, error_log, api_key):

    print 'Geocoding'

    # Read the backup csv file
    # This is used to make up for old places not being on wikidata
    institute_coordinates_backup = {}
    with open(config.config_dir + '/institute_coordinates.csv', 'rb') as csvfile:
        f = csv.reader(csvfile)
        for row in f:
            institute_coordinates_backup[row[0]] = {}
            institute_coordinates_backup[row[0]]['lat'] = row[1]
            institute_coordinates_backup[row[0]]['long'] = row[2]
    csvfile.close()

    locations_found = 0
    number_done = 1
    for this_paper in papers:
        found_coords = False

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

                # First we need to get the wikidata ID
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
                        this_paper['wikidata_item_id'] = item_id

                        # -- USE WIKIDATA ID TO GET GEO-COORDS
                        try:

                            retur = json.load(urllib2.urlopen('https://www.wikidata.org/w/api.php?action=wbgetentities&ids=' + item_id + '&format=json'))
                            p_lon = retur['entities'][item_id]['claims']['P625'][0]['mainsnak']['datavalue']['value']['longitude']
                            p_lat = retur['entities'][item_id]['claims']['P625'][0]['mainsnak']['datavalue']['value']['latitude']

                            locations_found += 1
                            this_paper['Extras']['LatLong'] = {'lat': str(p_lat), 'long': str(p_lon)}

                            found_coords = True

                        except:
                            try:

                                # There is not coordinate location property for the object Check headquaters location property
                                p_lon = retur['entities'][item_id]['claims']['P159'][0]['qualifiers']['P625'][0]['datavalue']['value']['longitude']
                                p_lat = retur['entities'][item_id]['claims']['P159'][0]['qualifiers']['P625'][0]['datavalue']['value']['latitude']

                                # print "Location Found from HQ " + papers[this_pmid]['Extras']['CleanInstitute']
                                locations_found += 1
                                this_paper['Extras']['LatLong'] = {'lat': str(p_lat), 'long': str(p_lon)}

                                found_coords = True

                            except:
                                print 'Unable to get geo-data (No HQ P625 Or P625) ' + this_paper['Extras']['CleanInstitute'] + " " + item_id + " (" + str(number_done) + "/" + str(len(papers)) + ")"

                    except:
                        # -- Use Backup File --
                        try:
                            backup_coords = institute_coordinates_backup[this_paper['Extras']['CleanInstitute']]
                            print 'Using Institute coordinates backup to get coordinates for ' + this_paper['Extras']['CleanInstitute']
                            this_paper['Extras']['LatLong'] = {'lat': str(backup_coords['lat']), 'long': str(backup_coords['long'])}
                            found_coords = True

                        except:
                            # The intitiue could not be found on wikidata or in the backup file
                            print 'Unable to get geo-data (Probably not on Wikidata) ' + this_paper['Extras']['CleanInstitute'] + " (" + str(number_done) + "/" + str(len(papers)) + ")"
                            error_log.logWarningPaper("Insititue " + this_paper['Extras']['CleanInstitute'] + " not on Wikidata (Consider using backup file)", this_paper)

                except:
                    print 'Unable to get geo-data (Wikidata Query Failed) ' + this_paper['Extras']['CleanInstitute'] + " (" + str(number_done) + "/" + str(len(papers)) + ")"
                    error_log.logWarningPaper("Wikidata query failed for " + this_paper['Extras']['CleanInstitute'], this_paper)

                # === End Look up ===

        except:
            print 'No Clean Institute for ' + this_paper['IDs']['hash'] + " (" + str(number_done) + "/" + str(len(papers)) + ")"
            error_log.logErrorPaper(" Clean Institute Missing for " + this_paper['IDs']['hash'] + " <a href='https://www.zotero.org/groups/300320/items/itemKey/" + this_paper['IDs']['zotero'] + "'>Zotero</a>", this_paper)

        # The coorditates have been found from either wikidata or the backup file
        # Now using the google maps api to get the country and city data
        if found_coords:
            # Get country for heatmap
            retur = json.load(urllib2.urlopen('https://maps.googleapis.com/maps/api/geocode/json?latlng=' + this_paper['Extras']['LatLong']['lat'] + ',' + this_paper['Extras']['LatLong']['long'] + '&key=' + api_key))

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
                error_log.logErrorPaper(this_paper['Extras']['CleanInstitute'] + " Google API Quota Reached!", this_paper)

        number_done += 1

    print "locations found: " + str(locations_found) + "/" + str(number_done)

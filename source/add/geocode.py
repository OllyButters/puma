#! /usr/bin/env python2

import csv
# import re
import os.path
import json
import urllib2
import requests
import logging

import config.config as config


# Have a go at geocoding the cleaned institute names
# I would expect there to be a lat long for all of them
def geocode(papers, error_log, api_key):

    print 'Geocoding'

    # Read in the backup csv file
    # This is used to make up for old places not being on wikidata or old places
    # being removed because it is difficult to reference them
    institute_coordinates_backup = {}
    with open(config.config_dir + '/institute_coordinates.csv', 'rb') as csvfile:
        f = csv.reader(csvfile)
        for row in f:
            institute_coordinates_backup[row[0]] = {}
            institute_coordinates_backup[row[0]]['lat'] = row[1]
            institute_coordinates_backup[row[0]]['long'] = row[2]
    csvfile.close()

    # Loop through the papers and attempt to find the coordinates, city name and country name
    locations_found = 0
    number_done = 0
    for this_paper in papers:

        found_coords = False

        try:

            # === Check if the location data is already cached ===
            # Each cached location is held in a seperate file
            clean = this_paper['clean']['location']['clean_institute']
            logging.info(str(this_paper['IDs']['hash']) + " geocode " + clean + " (" + str(number_done) + "/" + str(len(papers)) + ")")
            if os.path.isfile(config.cache_dir + "/geodata/" + clean):
                # The location has a cache file. load cached data.
                cache_file = open(config.cache_dir + "/geodata/" + clean, "r")
                data = cache_file.read()
                split = data.split("#")

                this_paper['clean']['location']['latitude'] = split[0]
                this_paper['clean']['location']['longitude'] = split[1]
                this_paper['clean']['location']['country_code'] = split[2]
                this_paper['clean']['location']['postal_town'] = split[3]

                locations_found += 1

            else:

                # === Look up clean institute coordinates on wikidata ===

                # First we need to get the wikidata ID
                # Form query to get the wikidata item by english label
                query = 'PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#> SELECT ?item WHERE { ?item rdfs:label "' + this_paper['clean']['location']['clean_institute'] + '"@en }'
                url = 'https://query.wikidata.org/bigdata/namespace/wdq/sparql'

                try:
                    data = requests.get(url, params={'query': query, 'format': 'json'}).json()

                    try:

                        # Process returned JSON to get entity id
                        item_uri = data['results']['bindings'][0]['item']['value']
                        item_uri_components = item_uri.split("/")
                        item_id = item_uri_components[len(item_uri_components)-1]
                        this_paper['wikidata_item_id'] = item_id

                        # === USE WIKIDATA ID TO GET GEO-COORDS ===
                        # The wikidata ID has been found now make a wikidata request for all the data
                        # about the item and look for the coordinate location
                        try:

                            retur = json.load(urllib2.urlopen('https://www.wikidata.org/w/api.php?action=wbgetentities&ids=' + item_id + '&format=json'))
                            p_lon = retur['entities'][item_id]['claims']['P625'][0]['mainsnak']['datavalue']['value']['longitude']
                            p_lat = retur['entities'][item_id]['claims']['P625'][0]['mainsnak']['datavalue']['value']['latitude']

                            locations_found += 1
                            this_paper['clean']['location']['latitude'] = str(p_lat)
                            this_paper['clean']['location']['longitude'] = str(p_lon)

                            found_coords = True

                        except:
                            try:
                                # The wikidata item does not have a statment for coordinate location. Instead check if there
                                # is coordinate location inside the headquaters location (if there is a HQ statment).
                                p_lon = retur['entities'][item_id]['claims']['P159'][0]['qualifiers']['P625'][0]['datavalue']['value']['longitude']
                                p_lat = retur['entities'][item_id]['claims']['P159'][0]['qualifiers']['P625'][0]['datavalue']['value']['latitude']

                                locations_found += 1
                                this_paper['clean']['location']['latitude'] = str(p_lat)
                                this_paper['clean']['location']['longitude'] = str(p_lon)

                                found_coords = True

                            except:
                                # No coordinate location or HQ coordinate location was found.
                                error_log.logWarningPaper('Unable to get geo-data (No HQ P625 Or P625) ' + this_paper['clean']['location']['clean_institute'] + " " + item_id + " (" + str(number_done) + "/" + str(len(papers)) + ")", this_paper)

                    except:
                        # === Not On Wikidata Check Backup File ===
                        # The item is not on wikidata. Check if the institute is in the backup coordinates file and get coordinates from there is possible.
                        try:
                            backup_coords = institute_coordinates_backup[this_paper['clean']['location']['clean_institute']]
                            # print 'Using Institute coordinates backup to get coordinates for ' + this_paper['Extras']['CleanInstitute']
                            this_paper['clean']['location']['latitude'] = str(backup_coords['lat'])
                            this_paper['clean']['location']['longitude'] = str(backup_coords['long'])

                            found_coords = True

                        except:
                            # The intitiue could not be found on wikidata or in the backup file
                            error_log.logWarningPaper("Insititue " + this_paper['clean']['location']['clean_institute'] + " not on Wikidata (Consider using backup file)", this_paper)

                except:
                    # Problem with the wikidata query
                    error_log.logWarningPaper("Wikidata query failed for " + this_paper['clean']['location']['clean_institute'], this_paper)

        except:
            # There is no clean institute for the paper object. Log this as an error so that it can be fixed in the Zotero file.
            logging.info('No Clean Institute for ' + this_paper['IDs']['hash'] + " (" + str(number_done) + "/" + str(len(papers)) + ")")
            error_log.logErrorPaper(" Clean Institute Missing for " + this_paper['IDs']['hash'] + " <a href='https://www.zotero.org/groups/" + config.zotero_id + "/items/itemKey/" + this_paper['IDs']['zotero'] + "'>Zotero</a>", this_paper)

        if found_coords:
            # The coordinates have been found from either wikidata or the backup file
            # Now using the google maps api to get the country and city data for use in the charts
            # Make API request.
            retur = json.load(urllib2.urlopen('https://maps.googleapis.com/maps/api/geocode/json?latlng=' + this_paper['clean']['location']['latitude'] + ',' + this_paper['clean']['location']['longitude'] + '&key=' + api_key))

            try:
                comps = retur['results'][0]['address_components']
                country_short = ""
                postal_town = ""
                for comp in comps:
                    # Search through the returned address components for the country and city names
                    if comp['types'][0] == "country":
                        country_short = comp['long_name']
                        this_paper['clean']['location']['country_code'] = country_short

                    if comp['types'][0] == "postal_town":
                        postal_town = comp['long_name']
                        this_paper['clean']['location']['postal_town'] = postal_town

                clean = clean.replace("/", "#")

                # Cache the data that has just been collected
                cache_file = open(config.cache_dir + "/geodata/" + clean, "w")
                cache_file.write(this_paper['clean']['location']['latitude'] + "#" + this_paper['clean']['location']['longitude'] + "#" + country_short + "#" + postal_town)
                cache_file.close()
            except:
                print 'Unable to get geo-data (Maybe Google API Quota Reached?) ' + this_paper['clean']['location']['clean_institute'] + " (" + str(number_done) + "/" + str(len(papers)) + ")"
                error_log.logErrorPaper(this_paper['clean']['location']['clean_institute'] + " Maybe Google API Quota Reached!", this_paper)

        number_done += 1

    print "locations found: " + str(locations_found) + "/" + str(number_done)

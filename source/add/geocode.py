#!/usr/bin/env python2

import csv
import os.path
import logging
import pprint

from SPARQLWrapper import SPARQLWrapper, JSON


import config.config as config


# Have a go at geocoding the cleaned institute names
# I would expect there to be a lat long for all of them
def geocode(papers, error_log, api_key):

    print('Geocoding')

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

    for this_paper in papers:

        found_coords = False

        # If there is no clean_institute then there really is no point trying this one.
        try:
            clean_institute = this_paper['clean']['location']['clean_institute']
        except:
            logging.warn('No clean_institute for ' + this_paper['IDs']['hash'])
            continue

        # Check if the location data is already cached.
        # Each cached location is held in a seperate file
        try:
            clean_institute = this_paper['clean']['location']['clean_institute']
            logging.info(str(this_paper['IDs']['hash']) + " geocode: " + clean_institute)
            if os.path.isfile(config.cache_dir + "/geodata/" + clean_institute):
                # The location has a cache file. load cached data.
                cache_file = open(config.cache_dir + "/geodata/" + clean_institute, "r")
                data = cache_file.read()
                split = data.split("#")

                this_paper['clean']['location']['latitude'] = split[0]
                this_paper['clean']['location']['longitude'] = split[1]
                this_paper['clean']['location']['country_code'] = split[2]
                this_paper['clean']['location']['postal_town'] = split[3]

                locations_found += 1
                logging.info("Added via cache file.")
                continue

        except:
            pass

        # Look up clean institute coordinates on wikidata
        try:
            # First we need to get the wikidata ID
            # Form query to get the wikidata item by english label
            # query = 'PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#> SELECT ?item WHERE { ?item rdfs:label "' + this_paper['clean']['location']['clean_institute'] + '"@en }'
            # url = 'https://query.wikidata.org/bigdata/namespace/wdq/sparql'
            url = 'https://query.wikidata.org/sparql'
            # v2
            # query = 'SELECT ?item WHERE { ?item rdfs:label "' + this_paper['clean']['location']['clean_institute'] + '"@en }'

            # v3

            query = '''
                SELECT ?item ?itemLabel ?country ?countryLabel ?mainTown ?mainTownLabel ?mainLon ?mainLat ?hqTownLabel ?hqLon ?hqLat
                WHERE {
                    ?item rdfs:label "''' + this_paper['clean']['location']['clean_institute'] + '''"@en.
                    ?item wdt:P17 ?country
                    OPTIONAL
                    {
                        ?item wdt:P131 ?mainTown
                    }
                    OPTIONAL
                    {
                        # Main location
                        ?item p:P625 ?mainLocation.
                        ?mainLocation psv:P625 ?mainCoordinateNode.
                        ?mainCoordinateNode wikibase:geoLongitude ?mainLon.
                        ?mainCoordinateNode wikibase:geoLatitude ?mainLat.
                    }
                    OPTIONAL
                    {
                        # HQ location
                        ?item wdt:P159 ?hqTown.
                    }
                    OPTIONAL
                    {
                        # HQ coordinates
                        ?item wdt:P159 ?hq2.
                        ?hq2 p:P625 ?hq3.
                        ?hq3 psv:P625 ?hq4.
                        ?hq4 wikibase:geoLongitude ?hqLon.
                        ?hq4 wikibase:geoLatitude ?hqLat.
                    }
                    SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
                }'''
            print(query)

            try:
                sparql = SPARQLWrapper(url)
                sparql.setQuery(query)
                sparql.setReturnFormat(JSON)
                data = sparql.query().convert()

                pp = pprint.PrettyPrinter(indent=4)

                pp.pprint(data)

                # print(data)
                logging.info(data)
                #
                # Parse the results
                # Country first
                try:
                    this_paper['clean']['location']['country'] = data['results']['bindings'][0]['countryLabel']['value']
                except:
                    pass

                # Try the main location
                try:
                    this_paper['clean']['location']['postal_town'] = data['results']['bindings'][0]['mainTownLabel']['value']
                    this_paper['clean']['location']['latitude'] = data['results']['bindings'][0]['mainLat']['value']
                    this_paper['clean']['location']['longitude'] = data['results']['bindings'][0]['mainLon']['value']
                except:
                    # Fall back to the HQ location if there is one
                    try:
                        this_paper['clean']['location']['postal_town'] = data['results']['bindings'][0]['hqTownLabel']['value']
                        this_paper['clean']['location']['latitude'] = data['results']['bindings'][0]['hqLat']['value']
                        this_paper['clean']['location']['longitude'] = data['results']['bindings'][0]['hqLon']['value']
                    except:
                        print("No suitable coordinates found from wikidata.")

                try:
                    print(this_paper['clean']['location']['country'])
                    print(this_paper['clean']['location']['postal_town'])
                    print(this_paper['clean']['location']['latitude'])
                    print(this_paper['clean']['location']['longitude'])
                except:
                    pass

            except Exception as e:
                # Problem with the wikidata query
                error_log.logWarningPaper("Wikidata query failed for " + this_paper['clean']['location']['clean_institute'], this_paper)
                logging.error("Wikidata query failed for " + this_paper['clean']['location']['clean_institute'] + str(e))
                print("Error with wikidata query")
                print(e)
        except:
            pass

        # Not On Wikidata Check Backup File
        if not found_coords:
            # The item is not on wikidata. Check if the institute is in the backup coordinates file and get coordinates from there is possible.
            try:
                backup_coords = institute_coordinates_backup[this_paper['clean']['location']['clean_institute']]
                # print 'Using Institute coordinates backup to get coordinates for ' + this_paper['Extras']['CleanInstitute']
                this_paper['clean']['location']['latitude'] = str(backup_coords['lat'])
                this_paper['clean']['location']['longitude'] = str(backup_coords['long'])

                locations_found += 1
                found_coords = True

            except:
                error_log.logWarningPaper("Insititue " + this_paper['clean']['location']['clean_institute'] + " not in backup file)", this_paper)

        # if found_coords and ('country_code' not in this_paper['clean']['location']) and ('postal_town' not in this_paper['clean']['location']):
            # The coordinates have been found from either wikidata or the backup file
            # but we don't have a country_code or postal_town
            # Now using the google maps api to get the country and city data for use in the charts
            # Make API request.
            # try:
            #    retur = json.load(urllib2.urlopen('https://maps.googleapis.com/maps/api/geocode/json?latlng=' + this_paper['clean']['location']['latitude'] + ',' + this_paper['clean']['location']['longitude'] + '&key=' + api_key))
            # except:
            #    print('Unable to get geo-data from Google API. ' + this_paper['clean']['location']['clean_institute'])

            # try:
            #    comps = retur['results'][0]['address_components']
            #    country_short = ""
            #    postal_town = ""
            #    for comp in comps:
            #        # Search through the returned address components for the country and city names
            #        if comp['types'][0] == "country":
            #            country_short = comp['long_name']
            #            this_paper['clean']['location']['country_code'] = country_short

            #        if comp['types'][0] == "postal_town":
            #            postal_town = comp['long_name']
            #            this_paper['clean']['location']['postal_town'] = postal_town

                # clean = clean.replace("/", "#")

        # Cache the data that has just been collected
        try:
            if 'latitude' in this_paper['clean']['location'] and 'longitude' in this_paper['clean']['location'] and 'country' in this_paper['clean']['location'] and 'postal_town' in this_paper['clean']['location']:
                cache_file = open(config.cache_dir + "/geodata/" + clean_institute, "w")
                cache_file.write(this_paper['clean']['location']['latitude'] + "#" + this_paper['clean']['location']['longitude'] + "#" + this_paper['clean']['location']['country'] + "#" + this_paper['clean']['location']['postal_town'])
                cache_file.close()
        except:
            pass
            # except:
            #    print("Error parsing Google geodata json file." + str(retur))
            #    error_log.logErrorPaper(this_paper['clean']['location']['clean_institute'] + " Maybe Google API Quota Reached!", this_paper)

    print("locations found: " + str(locations_found) + "/" + str(len(papers)))

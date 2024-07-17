
import csv
import os.path
import logging

from SPARQLWrapper import SPARQLWrapper, JSON

import config.config as config


# Have a go at geocoding the cleaned institute names
# I would expect there to be a lat long for all of them
def geocode(papers):

    print('Geocoding')
    logging.info('Geocoding.')

    # Read in the backup csv file
    # This is used to make up for old places not being on wikidata or old places
    # being removed because it is difficult to reference them
    institute_coordinates_backup = {}
    with open(config.config_dir + '/institute_coordinates.csv', encoding='utf8') as csvfile:
        f = csv.reader(csvfile)
        for row in f:
            institute_coordinates_backup[row[0]] = {}
            institute_coordinates_backup[row[0]]['lat'] = row[1]
            institute_coordinates_backup[row[0]]['long'] = row[2]
            try:
                institute_coordinates_backup[row[0]]['country'] = row[3]
            except Exception:
                pass
    csvfile.close()

    # Loop through the papers and attempt to find the coordinates and country name
    locations_found = 0

    for this_paper in papers:

        found_coords = False

        # If there is no clean_institute then there really is no point trying this one.
        try:
            clean_institute = this_paper['clean']['location']['clean_institute']
            logging.info(clean_institute)
        except Exception:
            logging.warning('No clean_institute for %s', this_paper['IDs']['hash'])
            continue

        # Check if the location data is already cached.
        # Each cached location is held in a separate file
        try:
            clean_institute = this_paper['clean']['location']['clean_institute']
            logging.info("%s geocode: %s", str(this_paper['IDs']['hash']), clean_institute)
            if os.path.isfile(config.cache_dir + "/geodata/" + clean_institute):
                # The location has a cache file. load cached data.
                cache_file = open(config.cache_dir + "/geodata/" + clean_institute, "r", encoding='utf8')
                data = cache_file.read()
                split = data.split("#")

                this_paper['clean']['location']['latitude'] = split[0]
                this_paper['clean']['location']['longitude'] = split[1]
                this_paper['clean']['location']['country'] = split[2]

                locations_found += 1
                logging.info("Added via cache file.")
                continue
        except Exception:
            pass

        # Look up clean institute coordinates on wikidata
        try:
            url = 'https://query.wikidata.org/sparql'

            query = '''
                SELECT ?item ?itemLabel ?country ?countryLabel ?mainTown ?mainTownLabel ?mainLon ?mainLat ?hqTownLabel ?hqLon ?hqLat
                WHERE {
                    ?item rdfs:label "''' + this_paper['clean']['location']['clean_institute'] + '''"@en.
                    ?item wdt:P17 ?country
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
                        # HQ coordinates
                        ?item wdt:P159 ?hq2.
                        ?hq2 p:P625 ?hq3.
                        ?hq3 psv:P625 ?hq4.
                        ?hq4 wikibase:geoLongitude ?hqLon.
                        ?hq4 wikibase:geoLatitude ?hqLat.
                    }
                    SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
                }'''

            logging.debug(query)

            try:
                sparql = SPARQLWrapper(url, agent="PUMA/0.1 (https://github.com/OllyButters/puma) Python3/SPARQLWrapper")
                sparql.setQuery(query)
                sparql.setReturnFormat(JSON)
                data = sparql.query().convert()

                # print(data)
                logging.info(data)

                # Parse the results
                # Country first
                try:
                    this_paper['clean']['location']['country'] = data['results']['bindings'][0]['countryLabel']['value']
                except Exception:
                    pass

                # Try the main location
                try:
                    this_paper['clean']['location']['latitude'] = data['results']['bindings'][0]['mainLat']['value']
                    this_paper['clean']['location']['longitude'] = data['results']['bindings'][0]['mainLon']['value']
                except Exception:
                    # Fall back to the HQ location if there is one
                    try:
                        this_paper['clean']['location']['latitude'] = data['results']['bindings'][0]['hqLat']['value']
                        this_paper['clean']['location']['longitude'] = data['results']['bindings'][0]['hqLon']['value']
                    except Exception:
                        logging.info("No suitable coordinates found from wikidata. (%s)", clean_institute)
                        print("No suitable coordinates found from wikidata. (" + clean_institute + ")")

                try:
                    logging.info(this_paper['clean']['location']['country'])
                except Exception:
                    pass

                # If I can print them then they must exist
                try:
                    logging.info(this_paper['clean']['location']['latitude'])
                    logging.info(this_paper['clean']['location']['longitude'])
                    found_coords = True
                except Exception:
                    pass

            except Exception as e:
                # Problem with the wikidata query
                logging.error("Wikidata query failed for %s %s",  this_paper['clean']['location']['clean_institute'], str(e))
                print("Error with wikidata query")
                print(e)
        except Exception:
            pass

        # Cache the data that has just been collected
        try:
            if 'latitude' in this_paper['clean']['location'] and 'longitude' in this_paper['clean']['location'] and 'country' in this_paper['clean']['location']:
                cache_file = open(config.cache_dir + "/geodata/" + clean_institute, "w", encoding='utf8')
                cache_file.write(this_paper['clean']['location']['latitude'] + "#" + this_paper['clean']['location']['longitude'] + "#" + this_paper['clean']['location']['country'])
                cache_file.close()
        except Exception:
            logging.warning("Problem caching geocode data for: %s", clean_institute)
            print("Problem caching geocode data for: " + clean_institute)

        # Not On Wikidata Check Backup File
        if not found_coords:
            # The item is not on wikidata. Check if the institute is in the backup coordinates file and get coordinates from there is possible.
            try:
                backup_coords = institute_coordinates_backup[this_paper['clean']['location']['clean_institute']]
                this_paper['clean']['location']['latitude'] = str(backup_coords['lat'])
                this_paper['clean']['location']['longitude'] = str(backup_coords['long'])
                this_paper['clean']['location']['country'] = str(backup_coords['country'])

                print("Coordinates added from backup file.")
                logging.info("Coordinates added from backup file.")
                locations_found += 1
                found_coords = True

            except Exception:
                logging.warning("Not in backup file either.")
                print("Not in backup file either.")

    print("locations found: " + str(locations_found) + "/" + str(len(papers)))

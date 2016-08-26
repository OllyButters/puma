#! /usr/bin/env python

import json
import re
import logging
import os
# import hashlib

import config.config as config


# Copy all the raw data to the processed directory, this means we are only
# ever working on the processed stuff and we never touch the raw data. This
# makes it easier to rerun as we don't have to rebuild the raw cache each time.
def pre_clean(papers, error_log):
    print 'precleaning'

    for this_paper in papers:

        # print this_paper['title']
        # Add IDs section
        this_paper['IDs'] = {}

        # Make hash from title
        # hash = hashlib.md5(this_paper['title'].encode('ascii', 'ignore')).hexdigest()
        # this_paper['IDs']['hash'] = hash

        # Ugly hack. Needs to be done better
        this_paper['IDs']['hash'] = this_paper['filename']

        this_paper['IDs']['DOI'] = ''
        this_paper['IDs']['PMID'] = ''
        this_paper['IDs']['zotero'] = ''

        # Add an extras item that we add stuff to - clean_institution,
        # citations etc
        this_paper['Extras'] = {}

        # Delete the empty authors.
        # There must be a more elegant way than this.
        authors_to_keep = []
        for i in range(0, len(this_paper['author'])-1):
            # print this_paper['author'][i]['family']
            try:
                if this_paper['author'][i]['family'] != "":
                    authors_to_keep.append(this_paper['author'][i])
            except:
                pass
        this_paper['author'] = authors_to_keep

        # Try sticking in the DOI
        try:
            this_paper['IDs']['DOI'] = this_paper['DOI']
        except:
            pass

        # Try sticking in the pmid
        try:
            this_paper['IDs']['PMID'] = this_paper['pmid']
        except:
            pass

        # Try sticking in the zotero id
        try:
            this_paper['IDs']['zotero'] = this_paper['key']
        except:
            pass

        # Generate a clean date
        # There are a lot of different dates in the paper data object.
        # These need to be convered into 1 date field so that it is consistently accessible throughout.
        # Relying on just this clean date will potentially cause some data lose so in situations where
        # you need a particular data, e.g. online publication date not paper publish date, then you
        # should make sure you are using the correct one.

        # CleanDate format = {'day':'00','month':'00','year':'0000'}
        this_paper['Extras']['CleanDate'] = {}

        try:

            if str(this_paper['PubmedData']['History'][0]['Day']) == "" or str(this_paper['PubmedData']['History'][0]['Month']) == "" or str(this_paper['PubmedData']['History'][0]['Year']) == "":
                raise Exception('Invalid Date')

            this_paper['Extras']['CleanDate']['day'] = str(this_paper['PubmedData']['History'][0]['Day'])
            this_paper['Extras']['CleanDate']['month'] = str(this_paper['PubmedData']['History'][0]['Month'])
            this_paper['Extras']['CleanDate']['year'] = str(this_paper['PubmedData']['History'][0]['Year'])

        except:
            try:
                if str(this_paper['issued']['date-parts'][0][2]) == "" or str(this_paper['issued']['date-parts'][0][1]) == "" or str(this_paper['issued']['date-parts'][0][0]) == "":
                    raise Exception('Invalid Date')

                this_paper['Extras']['CleanDate']['day'] = str(this_paper['issued']['date-parts'][0][2])
                this_paper['Extras']['CleanDate']['month'] = str(this_paper['issued']['date-parts'][0][1])
                this_paper['Extras']['CleanDate']['year'] = str(this_paper['issued']['date-parts'][0][0])

            except:
                try:
                    if str(this_paper['created']['date-parts'][0][2]) == "" or str(this_paper['created']['date-parts'][0][1]) == "" or str(this_paper['created']['date-parts'][0][0]) == "":
                        raise Exception('Invalid Date')

                    this_paper['Extras']['CleanDate']['day'] = str(this_paper['created']['date-parts'][0][2])
                    this_paper['Extras']['CleanDate']['month'] = str(this_paper['created']['date-parts'][0][1])
                    this_paper['Extras']['CleanDate']['year'] = str(this_paper['created']['date-parts'][0][0])

                except:
                    error_log.logErrorPaper("Cannot Create Clean Date", this_paper)


# Have a go at tidying up the mess that is first author institution.
# Essentially go through each institution and see if it matches a patten
# in the institute_cleaning.csv file. If it does then replace it with a
# standard name.
def clean_institution(papers):
    import unicodecsv
    logging.info('Starting institute cleaning')

    # Read in config file
    pattern = []
    replacements = []
    with open(config.config_dir + '/institute_cleaning.csv', 'rb') as csvfile:
        f = unicodecsv.reader(csvfile,  encoding='utf-8')
        for row in f:
            try:
                # Check it is not a comment string first.
                if re.match('#', row[0]):
                    continue

                # Check for blank lines
                if row[0] == '':
                    continue

                # If there is a second element in this row then carry on
                pattern.append(row[0])
                replacements.append(row[1])
            except:
                pass

    logging.info('Config read in, starting processing')

    # Cycle through institute checking the whole substitution list.
    # Stop when the first one matches.
    number_not_matched = 0
    for this_paper in papers:

        try:
            # print '============='
            # print this_paper['author'][0]['affiliation'][0]['name']
            institute = this_paper['author'][0]['affiliation'][0]['name']

        except:
            logging.warn('Could not find an affiliation for %s', this_paper)
            continue

        for y in range(0, len(pattern)):
            # logging.debug('%s %s %s', institute, pattern[y], replacements[y])
            temp = re.search(pattern[y], institute, re.IGNORECASE)
            if temp > 0:
                logging.info(
                    'ID:%s. %s MATCHES %s REPLACEDBY %s',
                    this_paper, institute, pattern[y], replacements[y])
                this_paper['Extras']['CleanInstitute'] = replacements[y]
                break

            if y == len(pattern)-1:
                logging.info('No match for %s. ID:%s', institute, this_paper)
                logging.warn('No match for %s. ID:%s', institute, this_paper)
                number_not_matched += 1

    print 'Cleaning institutions'
    print str(len(papers)-number_not_matched)+'/'+str(len(papers))+' cleaned'

    return number_not_matched


# Go through the deltas directory and apply any changes that are needed
def do_deltas(papers):

    delta_dir = config.config_dir + '/deltas/'

    deltas = os.listdir(delta_dir)

    print deltas

    for this_delta in deltas:
        # delta_file = '8680184'

        delta_path = delta_dir+this_delta

        fo = open(delta_path, 'r')
        record = json.loads(fo.read())
        fo.close()

        try:
            papers[this_delta]['Year'] = record['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['Year']
        except:
            print 'FAIL'

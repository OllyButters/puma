#! /usr/bin/env python

import json
import re
import logging
import os
# import hashlib

import config.config as config


# Some data will be missing. To deal with this, missing data will be put into the
# Zotero notes field ('extra') in a formatted structure "<key>:<value>\n<key>:<value>\n". That data will then be parsed by
# this function and put into the paper object so that it can be accessed throughout the analysis and html functions.
def clean_notes(papers, error_log):
    for this_paper in papers:
        try:
            notes = this_paper['merged']['extra']
            this_paper['merged']['extra'] = {}
            # split the notes data into key/value pairs
            fields = notes.split("\n")
            for this_field in fields:
                components = this_field.split(":")
                this_paper['extra'][components[0].strip()] = components[1].strip()
        except:
            pass


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
        for this_author in this_paper['author']:
            try:
                if this_author['family'] != "":
                    authors_to_keep.append(this_author)
            except:
                pass

        # for i in range(0, len(this_paper['author'])-1):
            # print this_paper['author'][i]['family']
        #     try:
        #         if this_paper['author'][i]['family'] != "":
        #            authors_to_keep.append(this_paper['author'][i])
        #    except:
        #        pass

        this_paper['author'] = authors_to_keep

        # add in a cleaned journal title
        clean_journal(this_paper)

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

        # Generate a Clean Date
        # There are a lot of different dates in the paper data object.
        # These need to be converted into 1 date field so that it is consistently accessible throughout.
        # Relying on just this clean date will potentially cause some data lose so in situations where
        # you need a particular date, e.g. online publication date not paper publish date, then you
        # should make sure you are using the correct one. However, the CleanDate field gives you the
        # best chance of getting a relevant date that is at least correct for something.

        # CleanDate format = {'day':'00','month':'00','year':'0000'}
        this_paper['Extras']['CleanDate'] = {}

        # Try the different date fields. If we don't get a full day, month and year for the CleanDate
        # then try the next possible field. Finally if none of the fields work then try the Zotero notes field.
        try:
            # First check for Pubmed date
            if str(this_paper['PubmedData']['History'][0]['Day']) == "" or str(this_paper['PubmedData']['History'][0]['Month']) == "" or str(this_paper['PubmedData']['History'][0]['Year']) == "":
                raise Exception('Invalid Date')

            this_paper['Extras']['CleanDate']['day'] = str(this_paper['PubmedData']['History'][0]['Day'])
            this_paper['Extras']['CleanDate']['month'] = str(this_paper['PubmedData']['History'][0]['Month'])
            this_paper['Extras']['CleanDate']['year'] = str(this_paper['PubmedData']['History'][0]['Year'])

        except:
            try:
                # Check for an issue date
                if str(this_paper['issued']['date-parts'][0][2]) == "" or str(this_paper['issued']['date-parts'][0][1]) == "" or str(this_paper['issued']['date-parts'][0][0]) == "":
                    raise Exception('Invalid Date')

                this_paper['Extras']['CleanDate']['day'] = str(this_paper['issued']['date-parts'][0][2])
                this_paper['Extras']['CleanDate']['month'] = str(this_paper['issued']['date-parts'][0][1])
                this_paper['Extras']['CleanDate']['year'] = str(this_paper['issued']['date-parts'][0][0])

            except:
                try:
                    # Check for a created date
                    if str(this_paper['created']['date-parts'][0][2]) == "" or str(this_paper['created']['date-parts'][0][1]) == "" or str(this_paper['created']['date-parts'][0][0]) == "":
                        raise Exception('Invalid Date')

                    this_paper['Extras']['CleanDate']['day'] = str(this_paper['created']['date-parts'][0][2])
                    this_paper['Extras']['CleanDate']['month'] = str(this_paper['created']['date-parts'][0][1])
                    this_paper['Extras']['CleanDate']['year'] = str(this_paper['created']['date-parts'][0][0])

                except:
                    try:
                        # Zotero Notes overide date
                        date_parts = this_paper['notes']['date'].split("/")
                        this_paper['Extras']['CleanDate']['day'] = str(date_parts[0])
                        this_paper['Extras']['CleanDate']['month'] = str(date_parts[1])
                        this_paper['Extras']['CleanDate']['year'] = str(date_parts[2])
                    except:
                        try:
                            # zotero 'date' field (only contains numerical year, word month)
                            date_parts = this_paper['date'].split(" ")
                            this_paper['Extras']['CleanDate']['year'] = str(date_parts[-1])
                        except:
                            # A date has not been found. Put this in the error log.
                            error_log.logErrorPaper("Cannot Create Clean Date (Consider using Zotero notes)", this_paper)


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
        f = unicodecsv.reader(csvfile,  encoding='utf-8')  # Handle extra unicode characters
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

        hasAffiliation = True
        try:
            # print '============='
            # print this_paper['author'][0]['affiliation'][0]['name']
            # institute = this_paper['PubmedArticle'][0]['MedlineCitation']['AuthorList'][0]['AffiliationInfo'][0]['Affiliation']
            institute = this_paper['author'][0]['affiliation'][0]['name']
            # institute = this_paper['PubmedArticle'][0]['MedlineCitation']['Article']['AuthorList'][0]['AffiliationInfo'][0]['Affiliation']
        except:
            logging.warn('Could not find an affiliation for %s', this_paper)
            hasAffiliation = False

        if hasAffiliation:
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

        try:
            this_paper['Extras']['CleanInstitute']
        except:
            # Check for Zotero note institution
            try:
                this_paper['Extras']['CleanInstitute'] = this_paper['notes']['institution']
            except:
                pass

    print 'Cleaning institutions'
    print str(len(papers)-number_not_matched) + '/' + str(len(papers)) + ' cleaned'

    return number_not_matched


# Clean journal title for paper
# Journal title /should/ be in this_paper['MedlineCitation']['Article']['Journal']['ISOAbbreviation'] but is sometimes missing. Look elsewhere (this_paper['container-title'] in addition.
# return this_paper with this_paper['cleaned-journal'] set
def clean_journal(this_paper):
    if not('cleaned-journal' in this_paper.keys() and this_paper['cleaned-journal'] is not None):
        this_paper['cleaned-journal'] = None
        try:
            this_paper['cleaned-journal'] = this_paper['MedlineCitation']['Article']['Journal']['ISOAbbreviation']
        except:
            try:
                this_paper['cleaned-journal'] = this_paper['container-title']
            except:
                pass
    return this_paper


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

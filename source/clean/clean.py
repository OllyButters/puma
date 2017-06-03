#! /usr/bin/env python

import json
import re
import logging
import os
import hashlib

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
                this_paper['merged']['extra'][components[0].strip()] = components[1].strip()
        except:
            pass


# Copy all the raw data to the processed directory, this means we are only
# ever working on the processed stuff and we never touch the raw data. This
# makes it easier to rerun as we don't have to rebuild the raw cache each time.
def pre_clean(papers, error_log):
    print 'precleaning'

    for this_paper in papers:

        # print this_paper['merged']['title']

        # Add an extras item that we add stuff to - clean_institution,
        # citations etc
        this_paper['clean'] = {}

        #add clean title
        this_paper['clean']['title'] = this_paper['merged']['title']

        # clean up the authors and add to 'clean'
        clean_authors(this_paper)

        # clean mesh headings into 'clean'
        clean_mesh(this_paper)

        # clean keywords into 'clean'
        clean_keywords(this_paper)

        # add in a cleaned journal title
        clean_journal(this_paper)

        # add in cleaned abstract
        clean_abstract(this_paper)

        # get a cleaned date
        # then insert year_published into clean['year_published'] if present
        clean_date(this_paper, error_log)
        try:
            this_paper['clean']['year_published'] = this_paper['clean']['cleaned_date']['year']
        except:
            pass

# get the abstract from MedlineCitation if present and add to clean['abstract']
def clean_abstract(this_paper):
    try:
        # Get abstract text
        this_paper['clean']['abstract'] = str(this_paper['merged']['MedlineCitation']['Article']['Abstract']['AbstractText'])
    except:
        logging.warn('No abstract for ' + this_paper['IDs']['hash'])

def clean_date(this_paper, error_log):
    # Generate a Clean Date
    # There are a lot of different dates in the paper data object.
    # These need to be converted into 1 date field so that it is consistently accessible throughout.
    # Relying on just this clean date will potentially cause some data lose so in situations where
    # you need a particular date, e.g. online publication date not paper publish date, then you
    # should make sure you are using the correct one. However, the clean_date field gives you the
    # best chance of getting a relevant date that is at least correct for something.

    # clean_date format = {'day':'00','month':'00','year':'0000'}
    this_paper['clean']['clean_date'] = {}

    # Try the different date fields. If we don't get a full day, month and year for the clean_date
    # then try the next possible field. Finally if none of the fields work then try the Zotero notes field.
    try:
        # First check for Pubmed date
        if str(this_paper['merged']['PubmedData']['History'][0]['Day']) == "" or str(this_paper['merged']['PubmedData']['History'][0]['Month']) == "" or str(this_paper['merged']['PubmedData']['History'][0]['Year']) == "":
            raise Exception('Invalid Date')

        this_paper['clean']['clean_date']['day'] = str(this_paper['merged']['PubmedData']['History'][0]['Day'])
        this_paper['clean']['clean_date']['month'] = str(this_paper['merged']['PubmedData']['History'][0]['Month'])
        this_paper['clean']['clean_date']['year'] = str(this_paper['merged']['PubmedData']['History'][0]['Year'])

    except:
        try:
            # Check for an issue date
            if str(this_paper['merged']['issued']['date-parts'][0][2]) == "" or str(this_paper['merged']['issued']['date-parts'][0][1]) == "" or str(this_paper['merged']['issued']['date-parts'][0][0]) == "":
                raise Exception('Invalid Date')

            this_paper['clean']['clean_date']['day'] = str(this_paper['merged']['issued']['date-parts'][0][2])
            this_paper['clean']['clean_date']['month'] = str(this_paper['merged']['issued']['date-parts'][0][1])
            this_paper['clean']['clean_date']['year'] = str(this_paper['merged']['issued']['date-parts'][0][0])

        except:
            try:
                # Check for a created date
                if str(this_paper['merged']['created']['date-parts'][0][2]) == "" or str(this_paper['merged']['created']['date-parts'][0][1]) == "" or str(this_paper['merged']['created']['date-parts'][0][0]) == "":
                    raise Exception('Invalid Date')

                this_paper['clean']['clean_date']['day'] = str(this_paper['merged']['created']['date-parts'][0][2])
                this_paper['clean']['clean_date']['month'] = str(this_paper['merged']['created']['date-parts'][0][1])
                this_paper['clean']['clean_date']['year'] = str(this_paper['merged']['created']['date-parts'][0][0])

            except:
                try:
                    # Zotero Notes overide date
                    date_parts = this_paper['merged']['notes']['date'].split("/")
                    this_paper['clean']['clean_date']['day'] = str(date_parts[0])
                    this_paper['clean']['clean_date']['month'] = str(date_parts[1])
                    this_paper['clean']['clean_date']['year'] = str(date_parts[2])
                except:
                    try:
                        # zotero 'date' field (only contains numerical year, word month)
                        date_parts = this_paper['merged']['date'].split(" ")
                        this_paper['clean']['clean_date']['year'] = str(date_parts[-1])
                    except:
                        # A date has not been found. Put this in the error log.
                        error_log.logErrorPaper("Cannot Create Clean Date (Consider using Zotero notes)", this_paper)

# Clean author name
# param author dict doi-style author object (at least 'family' and 'given' as keys)
# Author's cleaned name is lastname (family) followed by first initial
# returns string of cleaned name
def clean_author_name(this_author):
    # Create a clean author field. This is the Surname followed by first initial
    cleaned_author_name = this_author['family'] + " " + this_author['given'][0]
    return cleaned_author_name
      
# Clean up the author list
# Copies cleaned entry into this_paper['clean']['full_author_list']
# returns list of cleaned authors
def clean_authors(this_paper):
    # First, delete the empty authors.
    # There must be a more elegant way than this.
    authors_to_keep = []
    for this_author in this_paper['merged']['author']:
        try:
            if this_author['family'] != "":
                authors_to_keep.append(this_author)
        except:
            pass

    this_paper['merged']['author'] = authors_to_keep

    # now go through authors and clean name then append to clean full_author_list
    authors = []
    # Some pmid files dont actually have an authorlist! e.g. 2587412
    # This probably needs to be resolved with pubmed!
    try:
        # generate the relevant structure in clean
        this_paper['clean']['full_author_list'] = []
        
        for this_author in this_paper['merged']['author']:
            # There are some entries in the author list that are not actually authors e.g. 21379325 last author
            # note that any authors with an empty 'family' key have been
            # removed by clean.clean()
            try:
                authors.append(this_author['family'])
                # Create a clean author field. This is the Surname followed by first initial.
                #clean_author_name = this_author['family'] + " " + this_author['given'][0]
                cleaned_author_name = clean_author_name(this_author)
                this_author.update({'clean': cleaned_author_name})
                this_paper['clean']['full_author_list'].append(this_author)
            except:
                pass

    except Exception as e:
        print str(e)
        logging.warn('No AuthorList for ' + this_paper['IDs']['hash'])

    # do we need to return anything here? currently returns list of authors, probably should just be 0 or error code
    return authors

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

        #add location key in clean if not present
        if 'location' not in this_paper['clean'].keys():
          this_paper['clean']['location'] = {}

        hasAffiliation = True
        try:
            # print '============='
            # print this_paper['merged']['author'][0]['affiliation'][0]['name']
            # institute = this_paper['merged']['PubmedArticle'][0]['MedlineCitation']['AuthorList'][0]['AffiliationInfo'][0]['Affiliation']
            institute = this_paper['merged']['author'][0]['affiliation'][0]['name']
            # institute = this_paper['merged']['PubmedArticle'][0]['MedlineCitation']['Article']['AuthorList'][0]['AffiliationInfo'][0]['Affiliation']
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
                    this_paper['clean']['location']['clean_institute'] = replacements[y]

                    break

                if y == len(pattern)-1:
                    logging.info('No match for %s. ID:%s', institute, this_paper)
                    logging.warn('No match for %s. ID:%s', institute, this_paper)
                    number_not_matched += 1

        try:
            this_paper['clean']['location']['clean_institute']
        except:
            # Check for Zotero note institution
            try:
                this_paper['clean']['location']['clean_institute'] = this_paper['merged']['extra']['institution']
            except:
                pass

    print 'Cleaning institutions'
    print str(len(papers)-number_not_matched) + '/' + str(len(papers)) + ' cleaned'

    return number_not_matched


# Clean journal title for paper
# Journal title /should/ be in this_paper['merged']['MedlineCitation']['Article']['Journal']['ISOAbbreviation'] but is sometimes missing. Look elsewhere (this_paper['merged']['container-title'] in addition.
# return this_paper with this_paper['clean']['journal'] set
def clean_journal(this_paper):
    if not('journal' in this_paper['clean'].keys() and this_paper['clean']['journal'] is not None):
        this_paper['clean']['journal'] = {
          'journal_name': '',
          'volume': '',
          'issue': ''
        }
        try:
            this_paper['clean']['journal']['journal_name'] = this_paper['merged']['MedlineCitation']['Article']['Journal']['ISOAbbreviation']
        except:
            try:
                this_paper['clean']['journal']['journal_name'] = this_paper['merged']['container-title']
            except:
                logging.warn('No clean journal name for ' + this_paper['IDs']['hash'])
                pass
    return this_paper

# clean keywords (tags) for this_paper
def clean_keywords(this_paper):
    if 'tags' in this_paper['merged']:
        try:
            this_paper['clean']['keywords']
        except KeyError:
            this_paper['clean']['keywords'] = {}

        try:
            this_paper['clean']['keywords']['other']
        except KeyError:
            this_paper['clean']['keywords']['other'] = []

        try:
            for this_tag in this_paper['merged']['tags']:
                this_paper['clean']['keywords']['other'].append(
                    this_tag['tag'],
                )
        except:
            pass

# clean mesh headings for this_paper
def clean_mesh(this_paper):
    if 'MedlineCitation' in this_paper['merged'].keys():
        try:
            this_paper['clean']['keywords']
        except KeyError:
            this_paper['clean']['keywords'] = {}

        if 'MeshHeadingList' in this_paper['merged']['MedlineCitation'].keys():
            try:
                this_paper['clean']['keywords']['mesh']
            except KeyError:
                this_paper['clean']['keywords']['mesh'] = []

            try:
                for this_mesh in this_paper['merged']['MedlineCitation']['MeshHeadingList']:
                    this_paper['clean']['keywords']['mesh'].append(
                        {
                        'term': this_mesh['DescriptorName'],
                        'major': this_mesh['MajorTopicYN']
                        }
                    )
            except Exception as e:
                print str(e)
                pass

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
            papers[this_delta]['clean']['year_published'] = record['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['Year']
        except:
            print 'FAIL'

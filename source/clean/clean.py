#! /usr/bin/env python2

import re
import logging
import config.config as config


# Some data will be missing. To deal with this, missing data will be put into the
# Zotero notes field ('extra') in a formatted structure "<key>:<value>\n<key>:<value>\n". That data will then be parsed by
# this function and put into the paper object so that it can be accessed throughout the analysis and html functions.
def clean_notes(papers, error_log):
    for this_paper in papers:
        try:
            # notes = this_paper['merged']['extra']
            # this_paper['merged']['extra'] = {}
            notes = this_paper['raw']['zotero_data']['extra']
            this_paper['raw']['zotero_data']['extra'] = {}
            # split the notes data into key/value pairs
            fields = notes.split("\n")
            for this_field in fields:
                components = this_field.split(":")
                # this_paper['merged']['extra'][components[0].strip()] = components[1].strip()
                this_paper['clean']['zotero_data']['extra'][components[0].strip()] = components[1].strip()
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

        # add clean title
        this_paper['clean']['title'] = this_paper['raw']['zotero_data']['title']

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
# pmid_data   - usually present
# doi_data    - ??
# scopus_data - ??
def clean_abstract(this_paper):
    try:
        # Get abstract text
        this_paper['clean']['abstract'] = str(this_paper['raw']['pmid_data']['MedlineCitation']['Article']['Abstract']['AbstractText'])
    except:
        logging.warn('No abstract for ' + this_paper['IDs']['hash'])


def clean_date(this_paper, error_log):
    # Generate a Clean Date
    # There are a lot of different dates in the paper data object.
    # These need to be converted into 1 date field so that it is consistently accessible throughout.
    # Relying on just this clean date will potentially cause some data loss so in situations where
    # you need a particular date, e.g. online publication date not paper publish date, then you
    # should make sure you are using the correct one. However, the clean_date field gives you the
    # best chance of getting a relevant date that is at least correct for something.

    # clean_date format = {'day':'00','month':'00','year':'0000'}
    this_paper['clean']['clean_date'] = {}

    # Try the different date fields. If we don't get a full day, month and year for the clean_date
    # then try the next possible field. If that's not going well try getting just the year.
    # Finally if none of the fields work then try the Zotero notes field.
    date_status = False

    # Full Pubmed date - probably the best one
    if date_status is False:
        try:
            if str(this_paper['raw']['pmid_data']['PubmedData']['History'][0]['Day']) != "" and str(this_paper['raw']['pmid_data']['PubmedData']['History'][0]['Month']) != "" and str(this_paper['raw']['pmid_data']['PubmedData']['History'][0]['Year']) != "":
                this_paper['clean']['clean_date']['day'] = str(this_paper['raw']['pmid_data']['PubmedData']['History'][0]['Day'])
                this_paper['clean']['clean_date']['month'] = str(this_paper['raw']['pmid_data']['PubmedData']['History'][0]['Month'])
                this_paper['clean']['clean_date']['year'] = str(this_paper['raw']['pmid_data']['PubmedData']['History'][0]['Year'])
                date_status = True
        except:
            pass

    # Full issue date
    if date_status is False:
        try:
            if str(this_paper['raw']['pmid_data']['issued']['date-parts'][0][2]) != "" and str(this_paper['raw']['pmid_data']['issued']['date-parts'][0][1]) != "" and str(this_paper['raw']['pmid_data']['issued']['date-parts'][0][0]) != "":
                this_paper['clean']['clean_date']['day'] = str(this_paper['raw']['pmid_data']['issued']['date-parts'][0][2])
                this_paper['clean']['clean_date']['month'] = str(this_paper['raw']['pmid_data']['issued']['date-parts'][0][1])
                this_paper['clean']['clean_date']['year'] = str(this_paper['raw']['pmid_data']['issued']['date-parts'][0][0])
                date_status = True
        except:
            pass

    # Full created date
    # if date_status is False:
        # try:
            # actually I am not sure this is a good idea, I think this might be the date the metadata was created...
            # if str(this_paper['merged']['created']['date-parts'][0][2]) == "" or str(this_paper['merged']['created']['date-parts'][0][1]) == "" or str(this_paper['merged']['created']['date-parts'][0][0]) == "":
            #    raise Exception('Invalid Date')

            # this_paper['clean']['clean_date']['day'] = str(this_paper['merged']['created']['date-parts'][0][2])
            # this_paper['clean']['clean_date']['month'] = str(this_paper['merged']['created']['date-parts'][0][1])
            # this_paper['clean']['clean_date']['year'] = str(this_paper['merged']['created']['date-parts'][0][0])
            # date_status = True
        # except:
        #    pass

    # PARTIAL Pubmed date
    if date_status is False:
        try:
            if str(this_paper['raw']['pmid_data']['PubmedData']['History'][0]['Year']) != "":
                this_paper['clean']['clean_date']['year'] = str(this_paper['raw']['pmid_data']['PubmedData']['History'][0]['Year'])
                date_status = True
        except:
            pass

    # PARTIAL issue date
    if date_status is False:
        try:
            if str(this_paper['raw']['pmid_data']['issued']['date-parts'][0][0]) != "":
                this_paper['clean']['clean_date']['year'] = str(this_paper['raw']['pmid_data']['issued']['date-parts'][0][0])
                date_status = True
        except:
            pass

    # Zotero Notes overide date
    if date_status is False:
        try:
            date_parts = this_paper['raw']['zotero_data']['extra']['date'].split("/")
            this_paper['clean']['clean_date']['day'] = str(date_parts[0])
            this_paper['clean']['clean_date']['month'] = str(date_parts[1])
            this_paper['clean']['clean_date']['year'] = str(date_parts[2])
            date_status = True
        except:
            pass

    # zotero 'date' field (only contains numerical year, word month)
    if date_status is False:
        try:
            date_parts = this_paper['raw']['zotero_data']['date'].split(" ")
            this_paper['clean']['clean_date']['year'] = str(date_parts[-1])
            date_status = True
        except:
            pass

    # A date has not been found. Put this in the error log.
    if date_status is False:
        error_log.logErrorPaper("Cannot Create Clean Date (Consider using Zotero notes)", this_paper)
        logging.warn("Cannot Create Clean Date (Consider using Zotero notes). Hash: " + str(this_paper['IDs']['hash']))


# Clean author name
# param author dict doi-style author object (at least 'family' and 'given' as keys)
# Author's cleaned name is lastname (family) followed by first initial
# returns string of cleaned name
def clean_author_name(this_author):
    # Create a clean author field. This is the Surname followed by first initial
    cleaned_author_name = this_author['LastName'] + " " + this_author['Initials']
    return cleaned_author_name


# Clean up the author list
# Copies cleaned entry into this_paper['clean']['full_author_list']
# returns list of cleaned authors
# pmid_data   - usually present
# doi_data    - ??
# scopus_data - ??
def clean_authors(this_paper):
    # generate the relevant structure in clean
    this_paper['clean']['full_author_list'] = []

    # PMID first
    # Get the list of authors
    raw_author_list = []
    try:
        raw_author_list = this_paper['raw']['pmid_data']['MedlineCitation']['Article']['AuthorList']
    except:
        pass

    # First, delete the empty authors.
    # There must be a more elegant way than this.
    authors_to_keep = []
    # for this_author in this_paper['merged']['author']:
    for this_author in raw_author_list:
        try:
            if this_author['LastName'] != "":
                authors_to_keep.append(this_author)
        except:
            pass

    # this_paper['merged']['author'] = authors_to_keep

    # now go through authors and clean name then append to clean full_author_list
    authors = []
    # Some pmid files dont actually have an authorlist e.g. 2587412
    try:
        # for this_author in this_paper['merged']['author']:
        for this_author in authors_to_keep:
            # There are some entries in the author list that are not actually authors e.g. 21379325 last author
            # note that any authors with an empty 'family' key have been
            # removed by clean.clean()
            try:
                authors.append(this_author['LastName'])
                # Create a clean author field. This is the Surname followed by first initial.
                # clean_author_name = this_author['family'] + " " + this_author['given'][0]
                cleaned_author_name = clean_author_name(this_author)
                this_author.update({'clean': cleaned_author_name})
                this_paper['clean']['full_author_list'].append(this_author)
            except:
                pass

    except Exception as e:
        print str(e)
        logging.warn('No AuthorList for ' + this_paper['IDs']['hash'])

    # Sanity check we have something here, if not then see if the extras from zotero has something
    if len(this_paper['clean']['full_author_list']) == 0:
        logging.info('No author found, going to try extra field.')
        try:
            # this_paper['clean']['full_author_list'].append({'clean': this_paper['merged']['extra']['clean_first_author']})
            this_paper['clean']['full_author_list'].append({'clean': this_paper['raw']['zotero_data']['extra']['clean_first_author']})
            # this_paper['clean']['full_author_list'][0] = {'clean': this_paper['merged']['extra']['clean_first_author']}
            logging.debug('Author added via extra field.')
        except:
            logging.warn('Tried adding author from zotero extra, but failed. ')


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

                # Retype both to unicode, this will parse any \u1234 bits to their
                # actual unicode.
                # Make patter lowercase so it matches better.
                pattern.append(unicode(row[0]).lower())
                replacements.append(unicode(row[1]))
            except:
                pass

    # Stick a copy of the parsed lookup into the log.
    for i in range(0, len(pattern)):
        logging.debug(str(type(pattern[i])) + pattern[i] + " --> " + replacements[i])

    logging.info('Config read in, starting processing')

    # Cycle through institute checking the whole substitution list.
    # Stop when the first one matches.
    number_not_matched = 0
    for this_paper in papers:

        # add location key in clean if not present
        if 'location' not in this_paper['clean'].keys():
            this_paper['clean']['location'] = {}

        hasAffiliation = True
        try:
            institute = this_paper['raw']['pmid_data']['author'][0]['AffiliationInfo']['Affiliation']
        except:
            logging.warn('Could not find an affiliation for %s', this_paper['IDs']['zotero'])
            hasAffiliation = False

        if hasAffiliation:
            for y in range(0, len(pattern)):
                # logging.debug('%s %s %s', institute, pattern[y], replacements[y])

                # Check pattern in institite. These are both unicode and lowercase
                temp = pattern[y] in unicode(institute).lower()
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
                this_paper['clean']['location']['clean_institute'] = this_paper['merged']['extra']['clean_institute']
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
            this_paper['clean']['journal']['journal_name'] = this_paper['raw']['pmid_data']['MedlineCitation']['Article']['Journal']['ISOAbbreviation']
        except:
            try:
                this_paper['clean']['journal']['journal_name'] = this_paper['raw']['doi_data']['container-title']
            except:
                logging.warn('No clean journal name for ' + this_paper['IDs']['hash'])
                pass
    return this_paper


# clean keywords (tags) for this_paper
def clean_keywords(this_paper):
    if 'tags' in this_paper['raw']:
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
    if 'MedlineCitation' in this_paper['raw']['pmid_data'].keys():
        try:
            this_paper['clean']['keywords']
        except KeyError:
            this_paper['clean']['keywords'] = {}

        if 'MeshHeadingList' in this_paper['raw']['pmid_data']['MedlineCitation'].keys():
            try:
                this_paper['clean']['keywords']['mesh']
            except KeyError:
                this_paper['clean']['keywords']['mesh'] = []

            try:
                for this_mesh in this_paper['raw']['pmid_data']['MedlineCitation']['MeshHeadingList']:
                    this_paper['clean']['keywords']['mesh'].append(
                        {
                            'term': this_mesh['DescriptorName'],
                            'major': this_mesh['MajorTopicYN']
                        }
                    )
            except Exception as e:
                print str(e)
                pass

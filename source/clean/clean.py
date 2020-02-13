#!/usr/bin/env python3

import re
import logging
import hashlib
import config.config as config
import csv


################################################################################
# Copy and format all the relevant raw data into the clean part of the data object.
# This gives us a standard structure to work from later.
################################################################################
def clean(papers):
    print('Cleaning')
    logging.info('Starting cleaning.')

    for this_paper in papers:

        # Add an clean dict that we add stuff to
        this_paper['clean'] = {}

        # Parse the zotero extras field into clean/zoter_data
        parse_zotero_extras(this_paper)

        # clean the title
        clean_title(this_paper)

        # make a hash of the title
        this_paper['IDs']['hash'] = hashlib.md5(this_paper['clean']['title'].encode('ascii', 'ignore')).hexdigest()

        # clean up the authors and add to 'clean'
        clean_author_list(this_paper)

        # Figure out the first author (from the list or zotero)
        clean_first_author(this_paper)

        # clean mesh headings into 'clean'
        clean_mesh(this_paper)

        # clean keywords into 'clean'
        clean_keywords(this_paper)

        # add in a cleaned journal title
        clean_journal(this_paper)

        # add in cleaned abstract
        clean_abstract(this_paper)

        # get a cleaned date
        clean_date(this_paper)

        # Get the scopus data from the cache
        clean_citations_scopus(this_paper)

    logging.info('Finished cleaning.')
################################################################################


################################################################################
# Some data will be missing. To deal with this, missing data will be put into the
# Zotero extra field in a formatted structure "<key>:<value>\n<key>:<value>\n". That data will then be parsed by
# this function and put into the paper object so that it can be accessed throughout the analysis and html functions.
################################################################################
def parse_zotero_extras(this_paper):
    try:
        # Get the data, will bail here if there is none
        extras = this_paper['raw']['zotero_data']['extra']

        # Make the placeholders where this will go
        this_paper['clean']['zotero_data'] = {}
        this_paper['clean']['zotero_data']['extra'] = {}

        # split the extras data into key/value pairs
        fields = extras.split("\n")
        for this_field in fields:
            components = this_field.split(":")
            this_paper['clean']['zotero_data']['extra'][components[0].strip()] = components[1].strip()
    except:
        pass
################################################################################


################################################################################
# get the title and add to clean/title
# pmid_data   - usually present
# doi_data    - usually present
# scopus_data - usually present
# zotero_data - usually present
################################################################################
def clean_title(this_paper):
    clean_title_status = False

    # pmid first
    if not clean_title_status:
        try:
            this_paper['clean']['title'] = this_paper['raw']['pmid_data']['MedlineCitation']['Article']['ArticleTitle']
            clean_title_status = True
            logging.debug('Title added via PMID.' + this_paper['clean']['title'])
        except:
            pass

    # doi second
    if not clean_title_status:
        try:
            this_paper['clean']['title'] = this_paper['raw']['doi_data']['title']
            clean_title_status = True
            logging.debug('Title added via DOI.' + this_paper['clean']['title'])
        except:
            pass

    # scopus third
    if not clean_title_status:
        try:
            this_paper['clean']['title'] = this_paper['raw']['scopus_data']['search-results']['entry'][0]['dc:title']
            clean_title_status = True
            logging.debug('Title added via Scopus.' + this_paper['clean']['title'])
        except:
            pass

    # zotero last
    if not clean_title_status:
        try:
            this_paper['clean']['title'] = this_paper['raw']['zotero_data']['title']
            clean_title_status = True
            logging.debug('Title added via zotero.' + this_paper['clean']['title'])
        except:
            logging.warn('No abstract for ' + this_paper['IDs']['hash'])
################################################################################


################################################################################
# get the abstract from MedlineCitation if present and add to clean['abstract']
# pmid_data   - usually present
# doi_data    - never present
# scopus_data - never present
################################################################################
def clean_abstract(this_paper):
    try:
        # Get abstract text
        this_paper['clean']['abstract'] = str(this_paper['raw']['pmid_data']['MedlineCitation']['Article']['Abstract']['AbstractText'])
    except:
        logging.warn('No abstract for ' + this_paper['IDs']['hash'])
################################################################################


################################################################################
def clean_date(this_paper):
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

    ############################################################################
    # Try the various PubMed dates. Return True if one works
    def _clean_date_pmid(this_paper):

        # Just the year from the Journal info. This is what the journal wants, and may
        # be significantly later than when it first appeared online. It is how it would
        # appear in a citation.
        try:
            if str(this_paper['raw']['pmid_data']['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['Year']) != "":
                this_paper['clean']['clean_date']['year'] = str(this_paper['raw']['pmid_data']['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['Year'])
                return True
        except:
            pass

        # Full Pubmed date - I think this is derived from various places
        try:
            if str(this_paper['raw']['pmid_data']['MedlineCitation']['Article']['ArticleDate'][0]['Day']) != "" and str(this_paper['raw']['pmid_data']['MedlineCitation']['Article']['ArticleDate'][0]['Month']) != "" and str(this_paper['raw']['pmid_data']['MedlineCitation']['Article']['ArticleDate'][0]['Year']) != "":
                this_paper['clean']['clean_date']['day'] = str(this_paper['raw']['pmid_data']['MedlineCitation']['Article']['ArticleDate'][0]['Day'])
                this_paper['clean']['clean_date']['month'] = str(this_paper['raw']['pmid_data']['MedlineCitation']['Article']['ArticleDate'][0]['Month'])
                this_paper['clean']['clean_date']['year'] = str(this_paper['raw']['pmid_data']['MedlineCitation']['Article']['ArticleDate'][0]['Year'])
                return True
        except:
            pass

        # Full history
        try:
            if str(this_paper['raw']['pmid_data']['PubmedData']['History'][0]['Day']) != "" and str(this_paper['raw']['pmid_data']['PubmedData']['History'][0]['Month']) != "" and str(this_paper['raw']['pmid_data']['PubmedData']['History'][0]['Year']) != "":
                this_paper['clean']['clean_date']['day'] = str(this_paper['raw']['pmid_data']['PubmedData']['History'][0]['Day'])
                this_paper['clean']['clean_date']['month'] = str(this_paper['raw']['pmid_data']['PubmedData']['History'][0]['Month'])
                this_paper['clean']['clean_date']['year'] = str(this_paper['raw']['pmid_data']['PubmedData']['History'][0]['Year'])
                return True
        except:
            pass

        # Full issue date
        try:
            if str(this_paper['raw']['pmid_data']['issued']['date-parts'][0][2]) != "" and str(this_paper['raw']['pmid_data']['issued']['date-parts'][0][1]) != "" and str(this_paper['raw']['pmid_data']['issued']['date-parts'][0][0]) != "":
                this_paper['clean']['clean_date']['day'] = str(this_paper['raw']['pmid_data']['issued']['date-parts'][0][2])
                this_paper['clean']['clean_date']['month'] = str(this_paper['raw']['pmid_data']['issued']['date-parts'][0][1])
                this_paper['clean']['clean_date']['year'] = str(this_paper['raw']['pmid_data']['issued']['date-parts'][0][0])
                return True
        except:
            pass

        # PARTIAL Pubmed date
        try:
            if str(this_paper['raw']['pmid_data']['PubmedData']['History'][0]['Year']) != "":
                this_paper['clean']['clean_date']['year'] = str(this_paper['raw']['pmid_data']['PubmedData']['History'][0]['Year'])
                return True
        except:
            pass

        # PARTIAL issue date
        try:
            if str(this_paper['raw']['pmid_data']['issued']['date-parts'][0][0]) != "":
                this_paper['clean']['clean_date']['year'] = str(this_paper['raw']['pmid_data']['issued']['date-parts'][0][0])
                return True
        except:
            pass

        # if all fails, return False
        return False
    ############################################################################

    ############################################################################
    # Try the various DOI dates. Return True if one works
    def _clean_date_doi(this_paper):

        # journal-issue year
        try:
            if str(this_paper['raw']['doi_data']['journal-issue']['published-print']['date-parts'][0][0]) != "":
                this_paper['clean']['clean_date']['year'] = str(this_paper['raw']['doi_data']['journal-issue']['published-print']['date-parts'][0][0])
                return True
        except:
            pass

        # Try issued/data-parts (year and month)
        try:
            if str(this_paper['raw']['doi_data']['issued']['date-parts'][0][0]) != "" and str(this_paper['raw']['doi_data']['issued']['date-parts'][0][0]) != "":
                this_paper['clean']['clean_date']['month'] = str(this_paper['raw']['doi_data']['issued']['date-parts'][0][1])
                this_paper['clean']['clean_date']['year'] = str(this_paper['raw']['doi_data']['issued']['date-parts'][0][0])
                return True
        except:
            pass
        # Try issued/data-parts (year and month)
        try:
            if str(this_paper['raw']['doi_data']['issued']['date-parts'][0][0]) != "":
                this_paper['clean']['clean_date']['year'] = str(this_paper['raw']['doi_data']['issued']['date-parts'][0][0])
                return True
        except:
            pass

        # if all fails, return False
        return False

    ############################################################################
    # Try the Scopus data
    def _clean_date_scopus(this_paper):

        # check if only one result, otherwise return False
        try:
            if this_paper['raw']['scopus_data']['search-results']['opensearch:totalResults'] != '1':
                return False
        except:
            pass

        # most likely location is the prism:coverDate field
        try:
            cover_date_str = this_paper['raw']['scopus_data']['search-results']['entry'][0]['prism:coverDate']
            if cover_date_str != '':
                cover_date = cover_date_str.split('-')
                if len(cover_date[0]) == 4:
                    this_paper['clean']['clean_date']['year'] = cover_date[0]
                    if len(cover_date) == 3:
                        # assume iso date format (yyyy-mm-dd)
                        this_paper['clean']['clean_date']['month'] = cover_date[1]
                        this_paper['clean']['clean_date']['day'] = cover_date[2]
                    return True
        except:
            pass

        # if all fails, return False
        return False

    ############################################################################
    # Parse the zotero dates
    def _clean_date_zotero(this_paper):
        # Zotero Notes overide date
        try:
            date_parts = this_paper['raw']['zotero_data']['extra']['date'].split("/")
            this_paper['clean']['clean_date']['day'] = str(date_parts[0])
            this_paper['clean']['clean_date']['month'] = str(date_parts[1])
            this_paper['clean']['clean_date']['year'] = str(date_parts[2])
            return True
        except:
            pass

        # zotero 'date' field (only contains numerical year, word month)
        try:
            date_parts = this_paper['raw']['zotero_data']['date'].split(" ")
            this_paper['clean']['clean_date']['year'] = str(date_parts[-1])
            return True
        except:
            pass

        # if all fails, return False
        return False
    ############################################################################

    status = _clean_date_pmid(this_paper)
    if not status:
        status = _clean_date_doi(this_paper)
    if not status:
        status = _clean_date_scopus(this_paper)
    if not status:
        status = _clean_date_zotero(this_paper)

    # Make a note of the year published. This should be complete!
    if status:
        try:
            this_paper['clean']['year_published'] = this_paper['clean']['clean_date']['year']
        except:
            logging.warn("Cannot Create Clean Date (Consider using Zotero notes). Hash: " + str(this_paper['IDs']['hash']))

################################################################################


################################################################################
# Clean up the author list
# Copies cleaned entry into:
#
# this_paper['clean']['full_author_list']
# affiliation,
# given,
# clean (<Surname> <First initial>)
# family
#
# pmid_data   - usually present
# doi_data    - often present
# scopus_data - often present - NOT USED HERE.
################################################################################
def clean_author_list(this_paper):
    # generate the relevant structure in clean
    this_paper['clean']['full_author_list'] = []

    ############################################################################
    # PMID data
    ############################################################################
    def _clean_author_list_pmid(this_paper):
        try:
            # Get the list of authors
            raw_author_list = this_paper['raw']['pmid_data']['MedlineCitation']['Article']['AuthorList']

            # First, delete the empty authors.
            # There must be a more elegant way than this.
            authors = []
            for this_author in raw_author_list:
                try:
                    if this_author['LastName'] != "":
                        authors.append(this_author)
                except:
                    pass

            # now go through authors and clean name then append to clean full_author_list
            # Some pmid files dont actually have an authorlist e.g. 2587412
            try:
                for this_author in authors:
                    # There are some entries in the author list that are not actually authors e.g. 21379325 last author
                    # note that any authors with an empty 'LastName' key have been removed
                    try:
                        this_family = this_author['LastName']
                    except:
                        this_family = ''
                    try:
                        this_given = this_author['ForeName']
                    except:
                        this_given = ''
                    try:
                        this_affiliation = this_author['AffiliationInfo'][0]['Affiliation']
                    except:
                        this_affiliation = ''

                    if (this_family != '') and (this_given != ''):
                        this_clean = this_family + ' ' + this_given[0]

                    this_clean_author = {'clean': this_clean,
                                         'family': this_family,
                                         'given': this_given,
                                         'affiliation': {'name': this_affiliation}}

                    this_paper['clean']['full_author_list'].append(this_clean_author)
            except:
                pass

            # Check to see if we have something, if so we can
            if len(this_paper['clean']['full_author_list']) > 0:
                return True

        except:
            logging.warn('No AuthorList for ' + this_paper['IDs']['hash'])
    ############################################################################

    ############################################################################
    # DOI data
    ############################################################################
    def _clean_author_list_doi(this_paper):
        try:
            # Get the list of authors
            raw_author_list = this_paper['raw']['doi_data']['author']

            # First, delete the empty authors.
            # There must be a more elegant way than this.
            authors = []
            for this_author in raw_author_list:
                try:
                    if this_author['family'] != "":
                        authors.append(this_author)
                except:
                    pass

            # now go through authors and clean name then append to clean full_author_list
            try:
                for this_author in authors:
                    # There are some entries in the author list that are not actually authors e.g. 21379325 last author
                    # note that any authors with an empty 'LastName' key have been removed
                    try:
                        this_family = this_author['family']
                    except:
                        this_family = ''
                    try:
                        this_given = this_author['given']
                    except:
                        this_given = ''
                    try:
                        this_affiliation = this_author['affiliation'][0]['name']
                    except:
                        this_affiliation = ''

                    if (this_family != '') and (this_given != ''):
                        this_clean = this_family + ' ' + this_given[0]

                    this_clean_author = {'clean': this_clean,
                                         'family': this_family,
                                         'given': this_given,
                                         'affiliation': {'name': this_affiliation}}

                    this_paper['clean']['full_author_list'].append(this_clean_author)
            except:
                pass

            # Check to see if we have something, if so we can
            if len(this_paper['clean']['full_author_list']) > 0:
                return True

        except:
            logging.warn('No AuthorList for ' + this_paper['IDs']['hash'])
    ############################################################################

    # Actually run some code
    status = _clean_author_list_pmid(this_paper)

    if status is not True:
        status = _clean_author_list_doi(this_paper)
################################################################################


################################################################################
# The clean first author is not necessarily the first one in the clean author list
# it could have come from the zotero field clean_first_author, but that doesnt
# make sense to put in the list as it implies an author list of one which is
# probably not true.
################################################################################
def clean_first_author(this_paper):
    this_paper['clean']['first_author'] = ''
    try:
        # If the full_author_list has been populated (from any source) it will be here.
        this_paper['clean']['first_author'] = this_paper['clean']['full_author_list'][0]['clean']
    except:
        try:
            # try the zotero extras
            this_paper['clean']['first_author'] = this_paper['clean']['zotero_data']['extra']['clean_first_author']
        except:
            logging.warn('Tried adding author from zotero extra, but failed. ')
################################################################################


################################################################################
# Have a go at tidying up the mess that is first author institution.
# Essentially go through each institution and see if it matches a patten
# in the institute_cleaning.csv file. If it does then replace it with a
# standard name.
################################################################################
def clean_institution(papers):
    logging.info('Starting institute cleaning')

    # Read in config file
    pattern = []
    replacements = []
    with open(config.config_dir + '/institute_cleaning.csv', 'r') as inst_file:
        inst_file_reader = csv.reader(inst_file)
        # f = reader(csvfile,  encoding='utf-8')
        for row in inst_file_reader:
            logging.debug(row)
            try:
                # Check it is not a comment string first.
                if re.match('#', row[0]):
                    continue

                # Check for blank lines
                if row[0] == '':
                    continue

                # Retype both to unicode, this will parse any \u1234 bits to their
                # actual unicode.
                # Make pattern lowercase so it matches better.
                # pattern.append(unicode(row[0]).lower())
                # replacements.append(unicode(row[1]))
                pattern.append(str(row[0]).lower())
                replacements.append(str(row[1]))
            except Exception as e:
                print(e)
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
        if 'location' not in list(this_paper['clean'].keys()):
            this_paper['clean']['location'] = {}

        hasAffiliation = False

        #####
        # PubMed and DOI from first
        if not hasAffiliation:
            try:
                candidate_institute = this_paper['clean']['full_author_list'][0]['affiliation']['name']
                if candidate_institute != '':
                    hasAffiliation = True
            except:
                logging.info('No affil from PubMed for %s', this_paper['IDs']['zotero'])

        #####
        # Try scopus
        if not hasAffiliation:
            try:
                candidate_institute = this_paper['raw']['scopus_data']['search-results']['entry'][0]['affiliation'][0]['affilname']
                if candidate_institute != '':
                    hasAffiliation = True
            except:
                logging.info('No affil from Scopus for %s', this_paper['IDs']['zotero'])

        #####
        # Try zotero
        # Doing zotero here means it will still get passed through the matching
        # below to make sure it is a real place.
        if not hasAffiliation:
            try:
                candidate_institute = this_paper['clean']['zotero_data']['extra']['clean_institute']
                hasAffiliation = True
            except:
                logging.info('No affil from zotero for %s', this_paper['IDs']['zotero'])

        if hasAffiliation and candidate_institute != '':
            # Let's keep our candidate institite
            this_paper['clean']['location']['candidate_institute'] = candidate_institute

            for y in range(0, len(pattern)):
                # logging.debug('%s %s %s', institute, pattern[y], replacements[y])

                # Check pattern in institite. These are both unicode and lowercase
                # temp = pattern[y] in unicode(candidate_institute).lower()
                temp = pattern[y] in str(candidate_institute).lower()
                if temp > 0:
                    logging.info(
                        'ID:%s. %s MATCHES %s REPLACEDBY %s',
                        this_paper['IDs']['hash'], candidate_institute, pattern[y], replacements[y])
                    this_paper['clean']['location']['clean_institute'] = replacements[y]
                    break

                if y == len(pattern)-1:
                    logging.info('No match for %s. ID:%s', candidate_institute, this_paper['IDs']['hash'])
                    logging.warn('No match for %s. ID:%s', candidate_institute, this_paper['IDs']['hash'])
                    number_not_matched += 1

    print('Cleaning institutions')
    print(str(len(papers)-number_not_matched) + '/' + str(len(papers)) + ' cleaned')
################################################################################


################################################################################
# Clean journal title for paper
# return this_paper with this_paper['clean']['journal'] set
################################################################################
def clean_journal(this_paper):
    if not('journal' in list(this_paper['clean'].keys()) and this_paper['clean']['journal'] is not None):
        this_paper['clean']['journal'] = {
          'journal_name': '',
          'volume': '',
          'issue': ''
        }

        # PMID first
        try:
            candidate_journal = this_paper['raw']['pmid_data']['MedlineCitation']['Article']['Journal']['ISOAbbreviation']
            if isinstance(candidate_journal, str):
                this_paper['clean']['journal']['journal_name'] = candidate_journal
            elif isinstance(candidate_journal, list):
                this_paper['clean']['journal']['journal_name'] = candidate_journal[0]

            # Lets get the volume and issue too.
            try:
                this_paper['clean']['journal']['volume'] = this_paper['raw']['pmid_data']['MedlineCitation']['Article']['Journal']['JournalIssue']['Volume']
            except:
                pass

            try:
                this_paper['clean']['journal']['issue'] = this_paper['raw']['pmid_data']['MedlineCitation']['Article']['Journal']['JournalIssue']['Issue']
            except:
                pass

        # Try DOI if PMID didnt work
        except:
            try:
                candidate_journal = this_paper['raw']['doi_data']['container-title']
                if isinstance(candidate_journal, str):
                    this_paper['clean']['journal']['journal_name'] = candidate_journal
                elif isinstance(candidate_journal, list):
                    this_paper['clean']['journal']['journal_name'] = candidate_journal[0]

                # Lets get the volume and issue too.
                try:
                    this_paper['clean']['journal']['volume'] = this_paper['raw']['doi_data']['volume']
                except:
                    pass

                try:
                    this_paper['clean']['journal']['issue'] = this_paper['raw']['doi_data']['issue']
                except:
                    pass

            except:
                logging.warn('No clean journal name for ' + this_paper['IDs']['hash'])
                pass
################################################################################


################################################################################
# clean keywords (tags) for this_paper
# pmid - Usually present
# doi - never present
# scopus - never present
################################################################################
def clean_keywords(this_paper):
    try:
        this_paper['clean']['keywords']
    except KeyError:
        this_paper['clean']['keywords'] = {}

    try:
        this_paper['clean']['keywords']['other']
    except KeyError:
        this_paper['clean']['keywords']['other'] = []

    try:
        for this_tag in this_paper['raw']['pmid_data']['MedlineCitation']['KeywordList']:
            this_paper['clean']['keywords']['other'].append(
                this_tag['tag'],
            )
    except:
        pass
################################################################################


################################################################################
# clean mesh headings for this_paper
################################################################################
def clean_mesh(this_paper):
    try:
        this_paper['clean']['keywords']
    except KeyError:
        this_paper['clean']['keywords'] = {}

    try:
        this_paper['clean']['keywords']['mesh']
    except KeyError:
        this_paper['clean']['keywords']['mesh'] = []

    # Might not be any pmid data at all.
    try:
        if 'MeshHeadingList' in list(this_paper['raw']['pmid_data']['MedlineCitation'].keys()):
            try:
                for this_mesh in this_paper['raw']['pmid_data']['MedlineCitation']['MeshHeadingList']:
                    this_paper['clean']['keywords']['mesh'].append(
                        {
                            'term': this_mesh['DescriptorName'],
                            'major': this_mesh['MajorTopicYN']
                        }
                        )
            except Exception as e:
                logging.error('MeSH error: ' + str(e))
                print('MeSH error: ' + str(e))
    except:
        pass


################################################################################
# Get the scopus citation data from the cached scopus files.
################################################################################
def clean_citations_scopus(this_paper):
    # Do i really need this?
    this_paper['clean']['citations'] = {}
    this_paper['clean']['citations']['scopus'] = {}
    this_paper['clean']['citations']['PMC'] = {}

    try:
        # sometimes this returns multiple entries e.g. 22935244
        if len(this_paper['raw']['scopus_data']['search-results']['entry']) > 1:
            logging.warn("Multiple different citaton counts found!")
            return False

        if len(this_paper['raw']['scopus_data']['search-results']['entry']) == 1:
            this_paper['clean']['citations']['scopus']['count'] = this_paper['raw']['scopus_data']['search-results']['entry'][0]['citedby-count']
            logging.info('Citation added.')

            # Lets grab the scopus ID while we are here
            this_paper['IDs']['scopus'] = this_paper['raw']['scopus_data']['search-results']['entry'][0]['eid']
            return True

        if len(this_paper['raw']['scopus_data']['search-results']['entry']) == 0:
            logging.warn("0 citatons counts found!")
            return False
    except:
        pass
################################################################################

#! /usr/bin/env python2

import csv
import config.config as config
import os
import shutil
from nltk.stem import WordNetLemmatizer

############################################################
# Have all the data now, so do something with it
############################################################


# Build a list of all journals and count frequencies of each.
# output this to a csv file so it can be analysed by someone else.
def journals(papers):
    print "\n###Journals###"

    num_papers = len(papers)

    journals = []
    for this_paper in papers:
        if this_paper['clean']['journal']['journal_name'] != '':
            journals.append(this_paper['clean']['journal']['journal_name'])

    print str(len(journals)) + '/' + str(num_papers)
    print str(len(set(journals))) + ' different journals'

    # calculate the frequency of each journal
    freq = dict((x, journals.count(x)) for x in set(journals))

    i = 0
    print 'Top 5'

    # print a list of sorted frequencies
    with open(config.data_dir + '/journals.csv', 'wb') as csvfile:
        journals_file = csv.writer(csvfile)
        for w in sorted(freq, key=freq.get, reverse=True):
            if i < 5:
                print w, freq[w]
                i = i + 1
            journals_file.writerow([w.encode('utf-8'), freq[w]])


############################################################
# Build a list of words used in 'item'
############################################################
def word_frequencies(papers, item):
    print("\n###" + item + "###")

    # initialize lemmatizer
    lemmatizer = WordNetLemmatizer()

    all_words = []
    all_words_by_year = {}
    data_from_count = 0
    # Go through all papers
    for this_paper in papers:
        try:
            # Make sure there is a year dict for this year
            this_year = this_paper['clean']['year']
        except:
            this_year = '0'

        try:
            all_words_by_year[this_year]
        except:
            print(this_year)
            all_words_by_year[this_year] = {}

        try:
            # Get item text - this will be a long string of words
            text = str(this_paper['clean'][item])

            # Remove punctuation and esacpe characters that will cause a problem
            text = text.lower()
            text = text.replace("u'", "")
            text = text.replace("'", "")
            text = text.replace("{", "")
            text = text.replace("}", "")
            text = text.replace(",", " ")
            text = text.replace(".", " ")
            text = text.replace(":", " ")
            text = text.replace(";", " ")
            text = text.replace("'", "\'")
            text = text.replace('"', ' ')
            text = text.replace('[', ' ')
            text = text.replace(']', ' ')

            # The MeSH terms have things like major, and minor in them.
            if item == 'keywords':
                text = text.replace('mesh', '')
                text = text.replace('major', '')
                text = text.replace('minor', '')
                text = text.replace('term', '')

            # Get rid of the names of the studies. This should be in a config file REALLY.
            text = text.replace('alspac', '')
            text = text.replace('avon longitudinal study of pregnancy and childbirth', '')
            text = text.replace('avon longitudinal study of pregnancy and childhood', '')
            text = text.replace('avon longitudinal study of parents and children', '')
            text = text.replace('1958 birth cohort', '')
            text = text.replace('ncds', '')
            text = text.replace('national child development survey', '')

            # Add item words into list of all words
            all_words_by_year[this_year].extend(text.split())
            all_words.extend(text.split())
            data_from_count += 1
        except:
            pass

    # Parse all_words through a lemmatizer. This is like finding the stem, but
    # should always return real words.
    lemmatized_all_words = []
    lemmatized_all_words_by_year = {}

    for this_word in all_words:
        lemmatized_all_words.append(lemmatizer.lemmatize(this_word))

    for this_year in all_words_by_year:
        for this_word in all_words_by_year[this_year]:
            lemmatized_all_words_by_year[this_year].append(lemmatizer.lemmatize(this_word))

    # Stick this in a log file to check the lemmatizing makes sense
    with open(config.data_dir + '/' + item + '_lemmatizing_log.csv', 'wb') as csvfile:
        output_file = csv.writer(csvfile)
        for this_word in all_words:
            output_file.writerow([this_word, lemmatizer.lemmatize(this_word)])

    # calculate the frequency of each word in item
    raw_freq = dict((x, all_words.count(x)) for x in set(all_words))
    lemmatized_freq = dict((x, lemmatized_all_words.count(x)) for x in set(lemmatized_all_words))

    lemmatized_freq_by_year = {}
    for this_year in all_words_by_year:
        print(this_year)
        print(type(this_year))
        print(type(lemmatized_all_words_by_year))
        print(lemmatized_all_words_by_year)
        lemmatized_all_words_by_year[this_year] = {}
        print(lemmatized_all_words_by_year)
        bob = set(lemmatized_all_words_by_year[this_year])
        temp = dict((x, lemmatized_all_words_by_year[this_year].count(x)) for x in set(lemmatized_all_words_by_year[this_year]))
        # test = lemmatized_all_words_by_year[this_year]
        lemmatized_freq_by_year[this_year] = dict((x, lemmatized_all_words_by_year[this_year].count(x)) for x in set(lemmatized_all_words_by_year[this_year]))

    # = Remove stop words from the list of all words =
    # Read stop words from file
    stop_lines = tuple(open(config.config_dir + "/stopwords", "r"))
    stop_words = []
    for line in stop_lines:
        split = line.split()
        if len(split) > 0 and split[0] != "|" and "|" not in split[0]:
            stop_words.append(split[0])

    # Remove stop words
    for stp in stop_words:
        raw_freq.pop(stp, None)
        lemmatized_freq.pop(stp, None)

    i = 0
    print 'Top 5'

    # Output the RAW data to file and print the top 5 to screen.
    with open(config.data_dir + '/' + item + '_raw.csv', 'wb') as csvfile:
        output_file = csv.writer(csvfile)
        for w in sorted(raw_freq, key=raw_freq.get, reverse=True):
            if i < 5:
                print w, raw_freq[w]
                i = i+1
            output_file.writerow([w.encode('utf-8'), raw_freq[w]])

    # Output the LEMMATIZED data to file and print the top 5 to screen.
    with open(config.data_dir + '/' + item + '_lemmatized.csv', 'wb') as csvfile:
        output_file = csv.writer(csvfile)
        for w in sorted(lemmatized_freq, key=lemmatized_freq.get, reverse=True):
            if i < 5:
                print w, lemmatized_freq[w]
                i = i+1
            output_file.writerow([w.encode('utf-8'), lemmatized_freq[w]])

    # Output the LEMMATIZED data by year.
    # How do I actually want to output this? year as cols, words as rows? Will need to zero
    # all the cells first for when the word was not used in that year at all. Might need to run
    # through all the words in order first and build a zero dict.

    import numpy as np
    import pandas as pd

    # Years as columns
    # Words as rows

    # Get the list of years. This may have gaps
    years = all_words_by_year.keys()

    # years = ['62', '59', '58', '63', '67', '64']
    # words = ['study', 'chair', 'bobby', 'davey', 'house', 'thiss']

    # Make a list of years based on min and max years in the years list. It could be
    # that some years have no data. SS
    all_years = list(range(int(min(years)), int(max(years))+1))

    # Force the years to be strings, otherwise it is hard to index the cells if they are ints
    all_years = list(map(str, all_years))

    print(all_words)
    print(all_years)

    # Make a zero filled array
    len_years = int(len(all_years))
    len_words = int(len(all_words))
    A = np.zeros(len_years * len_words, dtype=int).reshape(len_words, len_years)

    # Make the DataFrame of the zeroes
    df = pd.DataFrame(A, index=sorted(all_words), columns=sorted(all_years))

    # Update a value in the DataFrame
    # df.at['study', '58'] = 123
    # df.at['house', '63'] = 123

    # for this_year in all_words_by_year:
    #     for this_word in all_words_by_year[this_year]:
    #         df.at[this_word,this_year] = freq

    # Update the DataFrame with the actual lemmatized data
    for this_year in lemmatized_freq_by_year:
        for this_word in lemmatized_freq_by_year[this_year]:
            df.at[this_word, this_year] = lemmatized_freq_by_year[this_year][this_word]

    print(len(lemmatized_freq_by_year))
    print(lemmatized_freq_by_year)

    # Output to a csv file
    df.to_csv('df.csv', index=True, header=True, sep=' ')

    exit(1)

    return data_from_count


############################################################
# Try with the authors - these are in a nested dict
############################################################
def authors(papers):

    import hashlib

    authors = []
    for this_paper in papers:
        try:
            for this_author in this_paper['clean']['full_author_list']:
                authors.append(this_author['clean'])
        except:
            pass

    # print authors
    freq = dict((x, authors.count(x)) for x in set(authors))
    print "\n###Authors###"

    print str(len(set(authors))) + ' different authors'
    # print freq

    i = 0
    print 'Top 5'

    with open(config.data_dir + '/authors.csv', 'wb') as csvfile:
        authors_file = csv.writer(csvfile)
        for w in sorted(freq, key=freq.get, reverse=True):
            if i < 5:
                print w, freq[w]
                i = i + 1
        # Need to utf-8 encode
            authors_file.writerow([w.encode('utf-8'), freq[w]])

    # === Author Network Object ===
    author_network = {}
    author_network['authors'] = {}
    author_network['connections'] = {}
    for this_paper in papers:
        try:
            # CREATE NODES
            for this_author in this_paper['clean']['full_author_list']:
                # Create author hash
                hash_object = hashlib.sha256(this_author)
                author_hash = hash_object.hexdigest()

                # Store author details
                try:
                    author_network['authors'][author_hash]
                except:
                    author_network['authors'][author_hash] = {}
                    author_network['authors'][author_hash]['clean'] = this_author['clean']

                try:
                    author_network['authors'][author_hash]['num_papers'] += 1
                except:
                    author_network['authors'][author_hash]['num_papers'] = 1

                # CREATE EDGES
                for con_author in this_paper['clean']['full_author_list']:

                    if not this_author == con_author:
                        con_hash_object = hashlib.sha256(con_author)
                        con_author_hash = con_hash_object.hexdigest()

                        con_id = ""
                        if author_hash > con_author_hash:
                            con_id = author_hash + "" + con_author_hash
                        else:
                            con_id = con_author_hash + "" + author_hash

                        con_id_hash_object = hashlib.sha256(con_id)
                        con_hash = con_id_hash_object.hexdigest()

                        try:
                            author_network['connections'][con_hash]
                        except:
                            author_network['connections'][con_hash] = {}
                            author_network['connections'][con_hash]['authors'] = []
                            author_network['connections'][con_hash]['authors'].append({'author_hash': author_hash})
                            author_network['connections'][con_hash]['authors'].append({'author_hash': con_author_hash})

                        try:
                            author_network['connections'][con_hash]['num_connections'] += 1
                        except:
                            author_network['connections'][con_hash]['num_connections'] = 1

        except:
            pass

    print "\n###Author network###"
    print str(len(author_network['authors'])) + " Authors"
    print str(len(author_network['connections'])) + " Connections"
    return author_network


############################################################
# Try with the FIRST authors - these are in a nested dict
############################################################
def first_authors(papers):

    num_papers = len(papers)
    first_authors = []
    for this_paper in papers:
        try:

            first_author_name = this_paper['clean']['full_author_list'][0]['clean']
            # stick the first author cleaned name in clean['first_author']['name']
            # try:
            #     this_paper['clean']['first_author']
            # except KeyError:
            #    this_paper['clean']['first_author'] = {}

            # this_paper['clean']['first_author']['name'] = first_author_name
            first_authors.append(first_author_name)

        except:
            next

    freq = dict((x, first_authors.count(x)) for x in set(first_authors))
    print "\n###First authors###"

    print str(len(first_authors))+'/'+str(num_papers)
    print str(len(set(first_authors)))+' different first authors'

    i = 0
    print 'Top 5'

    with open(config.data_dir + '/first_authors.csv', 'wb') as csvfile:
        authors_file = csv.writer(csvfile)
        for w in sorted(freq, key=freq.get, reverse=True):
            if i < 5:
                print w, freq[w]
                i = i+1
            # Need to utf-8 encode
            authors_file.writerow([w.encode('utf-8'), freq[w]])


############################################################
# Try with the FIRST authors INSTITUTE- these are in a nested dict
############################################################
def inst(papers):
    num_papers = len(papers)

    first_authors_inst = []
    for this_paper in papers:
        try:
            first_authors_inst.append(this_paper['clean']['location']['clean_institute'])
        except:
            pass

    freq = dict((x, first_authors_inst.count(x)) for x in set(first_authors_inst))
    print "\n###First authors institute###"

    print str(len(first_authors_inst))+'/'+str(num_papers)
    print str(len(set(first_authors_inst)))+' different first author institutes'

    i = 0
    print 'Top 5'

    with open(config.data_dir + '/first_authors_inst.csv', 'wb') as csvfile:
        authors_file = csv.writer(csvfile)
        for w in sorted(freq, key=freq.get, reverse=True):
            if i < 5:
                print w, freq[w]
                i = i+1
            # Need to utf-8 encode
            authors_file.writerow([w.encode('utf-8'), freq[w]])


############################################################
# Try with the mesh headings - these are in a nested dict
############################################################
def mesh(papers):
    num_papers = len(papers)

    mesh = []
    coverage = 0
    for this_paper in papers:

        if 'keywords' in this_paper['clean']['keywords'].keys():
            if 'mesh' in this_paper['clean']['keywords'].keys():
                coverage = coverage + 1

                try:
                    for this_mesh in this_paper['clean']['keywords']['mesh']:
                        mesh.append(this_mesh['term'])
                except:
                    pass

    freq = dict((x, mesh.count(x)) for x in set(mesh))
    print "\n###Mesh###"

    print str(coverage)+'/'+str(num_papers)
    print str(len(set(mesh)))+' different mesh headings'

    i = 0
    print 'Top 5'

    with open(config.data_dir + '/mesh.csv', 'wb') as csvfile:
        mesh_file = csv.writer(csvfile)
        for w in sorted(freq, key=freq.get, reverse=True):
            if i < 5:
                print w, freq[w]
                i = i+1
            # Need to utf-8 encode
            mesh_file.writerow([w.encode('utf-8'), freq[w]])


################################################################################
# output a csv file with some info in
################################################################################
def output_csv(papers):

    print '\n###Outputting CSV file###'
    with open(config.data_dir + '/all.csv', 'wb') as csvfile:
        all_file = csv.writer(csvfile)
        for this_paper in papers:
            try:
                title = this_paper['clean']['title']
            except:
                title = '???'

            try:
                first_author = this_paper['clean']['full_author_list'][0]['family']
            except:
                first_author = '???'

            try:
                journal = this_paper['clean']['journal']['journal_name']
            except:
                journal = '???'

            try:
                citations = this_paper['clean']['Citations']
            except:
                citations = '???'

            try:
                all_file.writerow([title, first_author, journal, citations, this_paper])
            except:
                pass


################################################################################
# Build an HTML report of the status of the important fields. Useful for cleaning
################################################################################
def coverage_report(papers):

    print 'Building coverage report'

    cov_css = '''
                .missing_required {background-color: red;}
                .missing_good_to_have {background-color: orange;}
                tr:nth-child(even) {background-color: lightgray}
                th {background-color: #4CAF50; color: white;}
                td, th {padding: 0.2em;}
                a:visited {color: red;}
    '''

    cov_html = '<table class="tablesorter">'
    cov_html += '''<thead><tr>
                    <th>Hash &uarr;&darr;</th>
                    <th>Zotero &uarr;&darr;</th>
                    <th>DOI &uarr;&darr;</th>
                    <th>PMID &uarr;&darr;</th>
                    <th>PMID DOI<br/>Lookup</th>
                    <th>PMID title<br/>Lookup</th>
                    <th>Title<br/>&uarr;&darr;</th>
                    <th>Keywords<br/>&uarr;&darr;</th>
                    <th>Abstract<br/>&uarr;&darr;</th>
                    <th>First<br/>Author &uarr;&darr;</th>
                    <th>First<br/>Author<br/>affil &uarr;&darr;</th>
                    <th>Clean<br/>Inst.<br/>&uarr;&darr;</th>
                    <th>Geocoded<br/>&uarr;&darr;</th>
                    <th>Clean<br/>Date &uarr;&darr;</th>
                    <th>Journal<br/>&uarr;&darr;</th>
                    <th>Scopus<br/>Citations &uarr;&darr;</th>
                </tr></thead>
                <tbody>'''

    status = {}
    status['hash'] = 0
    status['zotero'] = 0
    status['doi'] = 0
    status['pmid'] = 0
    status['title'] = 0
    status['keywords'] = 0
    status['abstract'] = 0
    status['first_author'] = 0
    status['first_author_affiliation'] = 0
    status['clean_institution'] = 0
    status['geocoded'] = 0
    status['clean_date'] = 0
    status['journal'] = 0
    status['scopus'] = 0

    for this_paper in papers:
        cov_html += '\n<tr class="item">'

        # Filename hash - this has to be prsent!
        try:
            fn_hash = this_paper['IDs']['hash']
            status['hash'] = status['hash'] + 1
        except:
            fn_hash = '???'

        cov_html += '<td><a href="status/cleaned/' + fn_hash + '" target="_blank">' + fn_hash + '</a></td>'

        # Zotero ID - this has to be present!
        try:
            zotero = this_paper['IDs']['zotero']
            cov_html += '<td><a href = "http://www.zotero.org/groups/' + config.zotero_id + '/items/itemKey/' + zotero + '" target="_blank">' + zotero + '</a></td>'
            status['zotero'] = status['zotero'] + 1
        except:
            cov_html += '<td class="missing_required">???</td>'

        # DOI - Not required, but REALLY useful
        try:
            doi = this_paper['IDs']['DOI']
            if doi != '':
                cov_html += '<td><a href="https://doi.org/' + doi + '" target="_blank">' + doi + '</a></td>'
                status['doi'] = status['doi'] + 1
            else:
                raise Exception()
        except:
            cov_html += '<td class="missing_good_to_have">???</td>'

        # PMID - Not required, but REALLY useful
        try:
            pmid = this_paper['IDs']['PMID']
            if pmid != '':
                cov_html += '<td><a href="https://www.ncbi.nlm.nih.gov/pubmed/' + pmid + '" target="_blank">' + pmid + '</a></td>'
                status['pmid'] = status['pmid'] + 1
            else:
                raise Exception()
        except:
            cov_html += '<td class="missing_good_to_have">???</td>'

        # PMID lookup based on DOI
        try:
            doi = this_paper['IDs']['DOI']
            if doi != '':
                # Swap the / for the html encoded version
                doi.replace('/', '%2F')
                cov_html += '<td><a href="https://www.ncbi.nlm.nih.gov/pubmed/?term=' + doi + '" target="_blank">Do lookup</a></td>'
            else:
                raise Exception()
        except:
            # Not bothered if not there
            cov_html += '<td></td>'

        # PMID lookup based on title
        try:
            title = this_paper['clean']['title']
            if title != '':
                ascii_title = this_paper['clean']['title'].encode("ascii", "ignore")
                cov_html += '<td><a href="https://www.ncbi.nlm.nih.gov/pubmed/?term=' + ascii_title + '" target="_blank">Do lookup</a></td>'
            else:
                raise Exception()
        except:
            # Not bothered if not there
            cov_html += '<td></td>'

        # Paper title
        try:
            title = this_paper['clean']['title']
            if title != '':
                cov_html += '<td>OK</td>'
                status['title'] = status['title'] + 1
            else:
                raise Exception()
        except:
            cov_html += '<td class="missing_good_to_have">???</td>'

        # Keywords
        try:
            keywords = this_paper['clean']['keywords']
            if keywords != '':
                cov_html += '<td>OK</td>'
                status['keywords'] = status['keywords'] + 1
            else:
                raise Exception()
        except:
            cov_html += '<td class="missing_good_to_have">???</td>'

        # Abstract
        try:
            abstract = this_paper['clean']['abstract']
            if abstract != '':
                cov_html += '<td>OK</td>'
                status['abstract'] = status['abstract'] + 1
            else:
                raise Exception()
        except:
            cov_html += '<td class="missing_good_to_have">???</td>'

        # First author - Not required, but REALLY useful
        try:
            first_author = this_paper['clean']['full_author_list'][0]['clean']
            if first_author != '':
                cov_html += '<td>OK</td>'
                status['first_author'] = status['first_author'] + 1
            else:
                raise Exception()
        except:
            cov_html += '<td class="missing_good_to_have">???</td>'

        # First author affiliation - Not required, but REALLY useful
        try:
            first_author_affiliation = this_paper['clean']['full_author_list'][0]['affiliation'][0]['name']
            if first_author_affiliation != '':
                cov_html += '<td>OK</td>'
                status['first_author_affiliation'] = status['first_author_affiliation'] + 1
            else:
                raise Exception()
        except:
            cov_html += '<td class="missing_good_to_have">???</td>'

        # CLEAN first author affiliation - Not required, but REALLY useful
        try:
            clean_institution = this_paper['clean']['location']['clean_institute']
            if clean_institution != '':
                cov_html += '<td>OK</td>'
                status['clean_institution'] = status['clean_institution'] + 1
            else:
                raise Exception()
        except:
            # Check to see if we have a DIRTY first author affiliation, if we do then we really should have a clean one too.
            try:
                first_author_affiliation = this_paper['clean']['full_author_list'][0]['affiliation'][0]['name']
                if first_author_affiliation != '':
                    cov_html += '<td class="missing_required">???</td>'
                else:
                    raise Exception()
            except:
                cov_html += '<td class="missing_good_to_have">???</td>'

        # Geocoded
        try:
            latitude = this_paper['clean']['location']['latitude']
            if latitude != '':
                cov_html += '<td>OK</td>'
                status['geocoded'] = status['geocoded'] + 1
            else:
                raise Exception()
        except:
            # Check to see if we have a clean institute, if we do then we should have a geocode
            try:
                clean_institution = this_paper['clean']['location']['clean_institute']
                if clean_institution != '':
                    cov_html += '<td class="missing_required">???</td>'
                else:
                    raise Exception()
            except:
                cov_html += '<td class="missing_good_to_have">???</!d>'

        # Clean date
        try:
            clean_date = this_paper['clean']['clean_date']['year']
            if clean_date != '':
                cov_html += '<td>' + this_paper['clean']['clean_date']['year'] + '</td>'
                status['clean_date'] = status['clean_date'] + 1
            else:
                raise Exception()
        except:
            cov_html += '<td class="missing_required">???</!d>'

        # Journal
        try:
            journal = this_paper['clean']['journal']['journal_name']
            if journal != '':
                cov_html += '<td>OK</td>'
                status['journal'] = status['journal'] + 1
            else:
                raise Exception()
        except:
            cov_html += '<td class="missing_good_to_have">???</td>'

        # Scopus citations
        try:
            scopus_url = 'http://api.elsevier.com/content/search/scopus'

            # Start a blank string to add to so can gracefully bail if problems
            scopus_html = ''

            # These should exist, they might be empty though
            doi = this_paper['IDs']['DOI']
            pmid = this_paper['IDs']['PMID']

            try:
                # Might not exist
                scopus = str(this_paper['clean']['citations']['scopus']['count'])

                # If there is some scopus data then show it
                if scopus != '':
                    scopus_html += '<td>'
                    scopus_html += scopus + '&nbsp;'
                    status['scopus'] = status['scopus'] + 1
            except:
                scopus_html += '<td class="missing_good_to_have">'

            # Put links in
            if doi != '':
                scopus_html += '<a href="' + scopus_url + '?apiKey=' + config.scopus_api_key + '&query=DOI(' + doi + ')&httpAccept=application%2Fjson" target="_blank">DOI</a>'
            if pmid != '':
                # Stick in a space if needed
                if doi != '':
                    scopus_html += '&nbsp;'
                scopus_html += '<a href="' + scopus_url + '?apiKey=' + config.scopus_api_key + '&query=PMID(' + pmid + ')&httpAccept=application%2Fjson" target="_blank">PMID</a>'

            scopus_html += '</td>'
            cov_html += scopus_html
        except:
            cov_html += '<td class="missing_required">ERROR</td>'

        cov_html += '</tr>'
    cov_html += '</tbody></table>'

    # Build a status table
    status_table = '<table>'
    status_table += '''<tr>
                    <th></th>
                    <th>Hash</th>
                    <th>Zotero</th>
                    <th>DOI</th>
                    <th>PMID</th>
                    <th>Title</th>
                    <th>Keywords</th>
                    <th>Abstract</th>
                    <th>First<br/>Author</th>
                    <th>First<br/>Author<br/>affil</th>
                    <th>Clean<br/>Inst</th>
                    <th>Geocoded</th>
                    <th>Clean<br/>Date</th>
                    <th>Journal</th>
                    <th>Scopus<br/>Citations</th>
                </tr>'''

    number_of_papers = len(papers)

    # Actual values
    status_table += '<tr>'
    status_table += '<td>Number (out of ' + str(number_of_papers) + ')</td>'
    status_table += '<td>' + str(status['hash']) + '</td>'
    status_table += '<td>' + str(status['zotero']) + '</td>'
    status_table += '<td>' + str(status['doi']) + '</td>'
    status_table += '<td>' + str(status['pmid']) + '</td>'
    status_table += '<td>' + str(status['title']) + '</td>'
    status_table += '<td>' + str(status['keywords']) + '</td>'
    status_table += '<td>' + str(status['abstract']) + '</td>'
    status_table += '<td>' + str(status['first_author']) + '</td>'
    status_table += '<td>' + str(status['first_author_affiliation']) + '</td>'
    status_table += '<td>' + str(status['clean_institution']) + '</td>'
    status_table += '<td>' + str(status['geocoded']) + '</td>'
    status_table += '<td>' + str(status['clean_date']) + '</td>'
    status_table += '<td>' + str(status['journal']) + '</td>'
    status_table += '<td>' + str(status['scopus']) + '</td>'
    status_table += '</tr>'

    # Percentages
    status_table += '<tr>'
    status_table += '<td>Percentage</td>'
    status_table += '<td>' + str(int(round(100*status['hash']/number_of_papers))) + '</td>'
    status_table += '<td>' + str(int(round(100*status['zotero']/number_of_papers))) + '</td>'
    status_table += '<td>' + str(int(round(100*status['doi']/number_of_papers))) + '</td>'
    status_table += '<td>' + str(int(round(100*status['pmid']/number_of_papers))) + '</td>'
    status_table += '<td>' + str(int(round(100*status['title']/number_of_papers))) + '</td>'
    status_table += '<td>' + str(int(round(100*status['keywords']/number_of_papers))) + '</td>'
    status_table += '<td>' + str(int(round(100*status['abstract']/number_of_papers))) + '</td>'
    status_table += '<td>' + str(int(round(100*status['first_author']/number_of_papers))) + '</td>'
    status_table += '<td>' + str(int(round(100*status['first_author_affiliation']/number_of_papers))) + '</td>'
    status_table += '<td>' + str(int(round(100*status['clean_institution']/number_of_papers))) + '</td>'
    status_table += '<td>' + str(int(round(100*status['geocoded']/number_of_papers))) + '</td>'
    status_table += '<td>' + str(int(round(100*status['clean_date']/number_of_papers))) + '</td>'
    status_table += '<td>' + str(int(round(100*status['journal']/number_of_papers))) + '</td>'
    status_table += '<td>' + str(int(round(100*status['scopus']/number_of_papers))) + '</td>'
    status_table += '</tr></table>'

    # Title
    title = '<h1>' + config.project_details['short_name'] + '</h1>'

    scripts = '<script src="https://ajax.googleapis.com/ajax/libs/jquery/3.3.1/jquery.min.js"></script>'
    scripts += '<script type="text/javascript" src="jquery.tablesorter.js"></script>'
    scripts += "<script>$(function(){$('table').tablesorter({widgets        : ['zebra', 'columns'],usNumberFormat : false,sortReset      : true,sortRestart    : true});});</script>"

    output_text = '<html><head><style>' + cov_css + '</style>' + scripts + '</head><body>' + title + status_table + '\n<br/><br/>' + cov_html + '</body></html>'
    coverage_file = open(config.html_dir + '/coverage_report.html', 'w')
    print >> coverage_file, output_text

    # Copy the jquery.tablesorter.js file across
    shutil.copy(config.template_dir + '/jquery.tablesorter.js', config.html_dir + '/jquery.tablesorter.js')

    # put a copy of all the processed files in the web tree
    if os.path.exists(config.html_dir + '/status/cleaned'):
        shutil.rmtree(config.html_dir + '/status/cleaned')
    shutil.copytree(config.cache_dir + '/processed/cleaned', config.html_dir + '/status/cleaned')

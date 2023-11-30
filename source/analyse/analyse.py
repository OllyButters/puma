#!/usr/bin/env python3

import csv
import config.config as config
from nltk.stem import WordNetLemmatizer
import logging
import numpy as np
import pandas as pd
import nltk
import re

# check if the 'wordnet' corpora for nltk is installed
# if not, download it
try:
    nltk.data.find('corpora/wordnet')
except:
    nltk.download('wordnet')

# check if the 'omw-1.4' corpora for nltk is installed
# if not, download it
try:
    nltk.data.find('corpora/omw-1.4')
except:
    nltk.download('omw-1.4')

############################################################
# Have all the data now, so do something with it
############################################################


############################################################
# Build a list of all journals and count frequencies of each.
# output this to a csv file so it can be analysed by someone else.
############################################################
def journals(papers):
    print("\n###Journals###")

    num_papers = len(papers)

    journals = []
    for this_paper in papers:
        if this_paper['clean']['journal']['journal_name'] != '':
            journals.append(this_paper['clean']['journal']['journal_name'])

    print(str(len(journals)) + '/' + str(num_papers))
    print(str(len(set(journals))) + ' different journals')

    # calculate the frequency of each journal
    freq = dict((x, journals.count(x)) for x in set(journals))

    i = 0
    print('Top 5')

    # print a list of sorted frequencies
    with open(config.data_dir + '/journals.csv', 'w') as csvfile:
        journals_file = csv.writer(csvfile)
        for w in sorted(freq, key=freq.get, reverse=True):
            if i < 5:
                print(w, freq[w])
                i = i + 1
            journals_file.writerow([w, freq[w]])


############################################################
# Build a list of words used in 'item'. Does this in a
# variety of ways:
# - Raw list of word frequencies
# - Lemmatized list of word frequencies
# - Lemmatized list of word frequencies by year
# - Lemmatized list of word frequencies by year weighted by
#   the number of papers published that year
############################################################
def word_frequencies(papers, item):
    print("\n### Doing word frequency on: " + item + "###")

    # Read stop words from file
    stop_lines = tuple(open(config.config_dir + "/stopwords", "r"))
    stop_words = []
    for line in stop_lines:
        stop_words.append(line.strip())

    logging.debug('Stop words:')
    logging.debug(stop_words)

    # initialize lemmatizer
    lemmatizer = WordNetLemmatizer()

    all_words = []    # A list
    all_words_by_year = {}  # A dictionary of lists
    number_of_papers_by_year = {}  # Needed to weight the values later
    data_from_count = 0
    # Go through all papers
    for this_paper in papers:
        try:
            # Make sure there is a year dict for this year
            this_year = this_paper['clean']['clean_date']['year']
        except:
            this_year = '0'

        # Increment the number of papers this year, or if we don't have one yet
        # initialize it to ONE!
        try:
            number_of_papers_by_year[this_year] = number_of_papers_by_year[this_year] + 1
        except:
            number_of_papers_by_year[this_year] = 1

        # Make sure there is a LIST for this year
        try:
            all_words_by_year[this_year]
        except:
            all_words_by_year[this_year] = []

        try:
            try:
                # Get item text - this will be a long string of words
                text = str(this_paper['clean'][item])
            except:
                logging.debug('No ' + item + ' text for this paper.')
                continue

            # Remove punctuation and escape characters that will cause a problem
            text = text.lower()
            text = text.replace("u'", '')
            text = text.replace('{', '')
            text = text.replace('}', '')
            text = text.replace('[', '')
            text = text.replace(']', '')
            text = text.replace('(', '')
            text = text.replace(')', '')
            text = text.replace(',', '')
            text = text.replace('.', '')
            text = text.replace(':', '')
            text = text.replace(';', '')
            text = text.replace("'", '')
            text = text.replace('"', '')
            text = text.replace('?', '')
            text = text.replace('!', '')

            # The MeSH terms have things like major, and minor in them.
            if item == 'keywords':
                text = text.replace('mesh', '')
                text = text.replace('major', '')
                text = text.replace('minor', '')
                text = text.replace('term', '')

            # Get rid of the names of the studies. This is done at this stage
            # as the whole string needs to be deleted, not individual words in itself.
            # This should be in a config file.
            text = text.replace('alspac', '')
            text = text.replace('avon longitudinal study of pregnancy and childbirth', '')
            text = text.replace('avon longitudinal study of pregnancy and childhood', '')
            text = text.replace('avon longitudinal study of parents and children', '')
            text = text.replace('children of the nineties', '')
            text = text.replace('1958 birth cohort', '')
            text = text.replace('ncds', '')
            text = text.replace('national child development survey', '')

            # split the text up for this paper using spaces
            temp_words = text.split()
            logging.debug('\n')
            logging.debug("Initial words: " + str(temp_words))

            # Only keep the word if it is just a word (i.e. not numbers or symbols)
            # and its not in the stop words file
            keep_words = []
            for this_word in temp_words:
                # Do a regex to make sure it uses a-z and optionally '-0-9
                if re.search("^[a-z]+['\-0-9]*$", this_word):
                    # Note that the output of lemmatizing may make a new word
                    # that is in the stopwords, but as that is run later
                    # it will not be removed. i.e. if you think this isn't
                    # working it is likely because the word is being added
                    # later on.
                    if this_word not in stop_words:
                        keep_words.append(this_word)
                    else:
                        logging.debug("Removing: " + str(this_word))
                else:
                    logging.debug("Removing: " + str(this_word))

            logging.debug("Keeping: " + str(keep_words))
            all_words.extend(keep_words)
            all_words_by_year[this_year].extend(keep_words)
            data_from_count += 1
        except Exception as e:
            logging.error('Uncaught error in word_frequencies: ' + str(e))
            print('Uncaught error in word_frequencies: ' + str(e))
            pass

    print('Total number of ' + item + ' words to keep: ' + str(len(all_words)))
    logging.debug('Total number of ' + item + ' words to keep: ' + str(len(all_words)))

    # Parse all_words through a lemmatizer. This is like finding the stem, but
    # should always return real words.
    lemmatized_all_words = []
    lemmatized_all_words_by_year = {}

    for this_word in all_words:
        lemmatized_all_words.append(lemmatizer.lemmatize(this_word))

    for this_year in all_words_by_year:
        # Make sure there is a LIST for this year
        try:
            lemmatized_all_words_by_year[this_year]
        except:
            lemmatized_all_words_by_year[this_year] = []

        # Do the lemmatization
        for this_word in all_words_by_year[this_year]:
            lemmatized_all_words_by_year[this_year].append(lemmatizer.lemmatize(this_word))

    # Stick this in a log file to check the lemmatizing makes sense
    with open(config.data_dir + '/' + item + '_lemmatizing_log.csv', 'w') as csvfile:
        output_file = csv.writer(csvfile)
        for this_word in all_words:
            output_file.writerow([this_word, lemmatizer.lemmatize(this_word)])

    # calculate the frequency of each word in item
    raw_freq = dict((x, all_words.count(x)) for x in set(all_words))
    lemmatized_freq = dict((x, lemmatized_all_words.count(x)) for x in set(lemmatized_all_words))

    lemmatized_freq_by_year = {}
    for this_year in all_words_by_year:
        lemmatized_freq_by_year[this_year] = dict((x, lemmatized_all_words_by_year[this_year].count(x)) for x in set(lemmatized_all_words_by_year[this_year]))

    # Output the RAW data to file.
    with open(config.data_dir + '/' + item + '_raw.csv', 'w') as csvfile:
        output_file = csv.writer(csvfile)
        for w in sorted(raw_freq, key=raw_freq.get, reverse=True):
            output_file.writerow([w, raw_freq[w]])

    i = 0
    print('Top 5 Lemmatized words:')

    # Output the LEMMATIZED data to file and print the top 5 to screen.
    with open(config.data_dir + '/' + item + '_lemmatized.csv', 'w') as csvfile:
        output_file = csv.writer(csvfile)
        for w in sorted(lemmatized_freq, key=lemmatized_freq.get, reverse=True):
            if i < 5:
                print(w, lemmatized_freq[w])
                i = i+1
            output_file.writerow([w, lemmatized_freq[w]])

    # Output the LEMMATIZED data by year.
    # Years as columns
    # Words as rows

    # Get the list of years. This may have gaps
    years = list(all_words_by_year.keys())

    # Might be a zero year if a year was not set in the clean data.
    # Drop the zero for now, otherwise we have a list from 0-max(years).
    try:
        years.remove('0')
    except:
        pass

    # Make a list of years based on min and max years in the years list. It could be
    # that some years have no data.
    print('Min/max years: ' + str(int(min(years))) + '/' + str(int(max(years))))
    all_years = list(range(int(min(years)), int(max(years))+1))

    # Force the years to be strings, otherwise it is hard to index the cells if they are ints
    all_years = list(map(str, all_years))

    # Make a zero filled array
    len_years = int(len(set(all_years)))
    len_words = int(len(set(lemmatized_all_words)))
    print('Unique Lemmatized words: ' + str(len_words))
    A = np.zeros(len_years * len_words, dtype=int).reshape(len_words, len_years)

    # Make the DataFrame of the zeroes
    df = pd.DataFrame(A, index=sorted(set(lemmatized_all_words)), columns=sorted(set(all_years)))

    # Update the DataFrame with the actual lemmatized data
    for this_year in lemmatized_freq_by_year:
        for this_word in lemmatized_freq_by_year[this_year]:
            df.at[this_word, this_year] = lemmatized_freq_by_year[this_year][this_word]

    # Output to a csv file
    df.to_csv(config.data_dir + '/' + item + '_lemmatized_by_year.csv', index=True, header=True, sep=',')

    # Do the same but weight it by the number of papers published this year
    # Make the DataFrame of the zeroes
    A = np.zeros(len_years * len_words, dtype=float).reshape(len_words, len_years)
    df_weighted = pd.DataFrame(A, index=sorted(set(lemmatized_all_words)), columns=sorted(set(all_years)))

    # Update the DataFrame with the actual lemmatized data
    for this_year in lemmatized_freq_by_year:
        for this_word in lemmatized_freq_by_year[this_year]:
            # Put in zero check to stop divide by zero errors. If the number_of_papers_by_year is zero then
            # the lemmatized_freq_by_year for this year should also be zero, so this checks if that is not the case.
            if number_of_papers_by_year[this_year] == 0:
                if lemmatized_freq_by_year[this_year][this_word] == 0:
                    continue
                else:
                    print('ERROR! No papers published this year, but there are words for this year! ' + this_year)
                    print(number_of_papers_by_year)
                    print(lemmatized_freq_by_year[this_year])
                    exit(1)
            df_weighted.at[this_word, this_year] = (1.0*lemmatized_freq_by_year[this_year][this_word])/(1.0*number_of_papers_by_year[this_year])

    # Output to a csv file
    df_weighted.to_csv(config.data_dir + '/' + item + '_lemmatized_by_year_weighted.csv', index=True, header=True, sep=',')

    return data_from_count


############################################################
# Calculate frequency of all authors and output to file
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
    print("\n###Authors###")

    print(str(len(set(authors))) + ' different authors')
    # print freq

    i = 0
    print('Top 5')

    with open(config.data_dir + '/authors.csv', 'w') as csvfile:
        authors_file = csv.writer(csvfile)
        for w in sorted(freq, key=freq.get, reverse=True):
            if i < 5:
                print(w, freq[w])
                i = i + 1
            authors_file.writerow([w, freq[w]])

    # === Author Network Object ===
    author_network = {}
    author_network['authors'] = {}
    author_network['connections'] = {}
    for this_paper in papers:
        try:
            # CREATE NODES
            for this_author in this_paper['clean']['full_author_list']:
                # Create author hash

                # hash_object = hashlib.sha256(this_author)
                # author_hash = hash_object.hexdigest()

                author_hash = hashlib.sha256(this_author['clean'].encode('ascii', 'ignore')).hexdigest()

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
                        # con_hash_object = hashlib.sha256(con_author)
                        # con_author_hash = con_hash_object.hexdigest()

                        con_author_hash = hashlib.sha256(con_author['clean'].encode('ascii', 'ignore')).hexdigest()

                        con_id = ""
                        if author_hash > con_author_hash:
                            con_id = author_hash + "" + con_author_hash
                        else:
                            con_id = con_author_hash + "" + author_hash

                        # con_id_hash_object = hashlib.sha256(con_id)
                        # con_hash = con_id_hash_object.hexdigest()
                        con_hash = hashlib.sha256(con_id.encode('ascii', 'ignore')).hexdigest()

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

        except Exception as e:
            print('Uncaught error: ' + str(e))
            pass

    print("\n###Author network###")
    print(str(len(author_network['authors'])) + " Authors")
    print(str(len(author_network['connections'])) + " Connections")
    return author_network


############################################################
# Calculate frequency of first authors and output to file
############################################################
def first_authors(papers):

    num_papers = len(papers)
    first_authors = []
    for this_paper in papers:
        try:
            first_author_name = this_paper['clean']['first_author']
            first_authors.append(first_author_name)
        except:
            next

    freq = dict((x, first_authors.count(x)) for x in set(first_authors))
    print("\n###First authors###")

    print(str(len(first_authors))+'/'+str(num_papers))
    print(str(len(set(first_authors)))+' different first authors')

    i = 0
    print('Top 5')

    with open(config.data_dir + '/first_authors.csv', 'w') as csvfile:
        authors_file = csv.writer(csvfile)
        for w in sorted(freq, key=freq.get, reverse=True):
            if i < 5:
                print(w, freq[w])
                i = i+1
            authors_file.writerow([w, freq[w]])


############################################################
# Calculate frequency of first authors' institutes and output to file
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
    print("\n###First authors institute###")

    print(str(len(first_authors_inst))+'/'+str(num_papers))
    print(str(len(set(first_authors_inst)))+' different first author institutes')

    i = 0
    print('Top 5')

    with open(config.data_dir + '/first_authors_inst.csv', 'w') as csvfile:
        authors_file = csv.writer(csvfile)
        for w in sorted(freq, key=freq.get, reverse=True):
            if i < 5:
                print(w, freq[w])
                i = i+1
            # Need to utf-8 encode
            authors_file.writerow([w, freq[w]])


############################################################
# Calculate frequency of all mesh and output to file
############################################################
def mesh(papers):
    num_papers = len(papers)

    mesh = []
    coverage = 0
    for this_paper in papers:

        if 'keywords' in list(this_paper['clean']['keywords'].keys()):
            if 'mesh' in list(this_paper['clean']['keywords'].keys()):
                coverage = coverage + 1

                try:
                    for this_mesh in this_paper['clean']['keywords']['mesh']:
                        mesh.append(this_mesh['term'])
                except:
                    pass

    freq = dict((x, mesh.count(x)) for x in set(mesh))
    print("\n###Mesh###")

    print(str(coverage)+'/'+str(num_papers))
    print(str(len(set(mesh)))+' different mesh headings')

    i = 0
    print('Top 5')

    with open(config.data_dir + '/mesh.csv', 'w') as csvfile:
        mesh_file = csv.writer(csvfile)
        for w in sorted(freq, key=freq.get, reverse=True):
            if i < 5:
                print(w, freq[w])
                i = i+1
            # Need to utf-8 encode
            mesh_file.writerow([w, freq[w]])


################################################################################
# output a csv file with some general info in
################################################################################
def output_csv(papers):

    print('\n###Outputting CSV file###')
    with open(config.data_dir + '/all.csv', 'w') as csvfile:
        all_file = csv.writer(csvfile)
        all_file.writerow(['year', 'authors', 'title', 'first_author', 'journal', 'citations', 'tags'])

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
                all_authors = this_paper['clean']['full_author_list']
                author_string = ''
                for this_author in all_authors:
                    author_string = author_string + this_author['family'] + ', '
            except:
                author_string = '???'

            try:
                journal = this_paper['clean']['journal']['journal_name']
            except:
                journal = '???'

            try:
                citations = this_paper['clean']['citations']['scopus']['count']
            except:
                citations = '???'

            try:
                year = this_paper['clean']['clean_date']['year']
            except:
                year = '???'

            try:
                tags = []
                tag_string = ''
                for this_tag in this_paper['clean']['zotero_tags']:
                    tags.append(this_tag['tag'])
                tag_string = ', '.join(tags)
            except:
                tag_string = '???'

            try:
                all_file.writerow([year, author_string, title, first_author, journal, citations, tag_string])
            except:
                pass

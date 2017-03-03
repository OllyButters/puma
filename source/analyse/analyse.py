#! /usr/bin/env python

import csv
import logging

import config.config as config

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
        try:
            journals.append(this_paper['journalAbbreviation'])
        except:
            logging.warn('No journalAbbreviation for '+this_paper['IDs']['hash'])

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
# Build a list of words used in Abstracts
############################################################
def abstracts(papers):
    print "\n###Abstracts###"

    words = []
    data_from_count = 0
    # Go through all papers
    for this_paper in papers:
        try:
            # Get abstract text
            abstracts = str(this_paper['MedlineCitation']['Article']['Abstract']['AbstractText'])

            # Remove punctuation and esacpe characters that will cause a problem
            abstracts = abstracts.lower()
            abstracts = abstracts.replace(",", " ")
            abstracts = abstracts.replace(".", " ")
            abstracts = abstracts.replace(":", " ")
            abstracts = abstracts.replace(";", " ")
            abstracts = abstracts.replace("'", "\'")
            abstracts = abstracts.replace('"', '\"')

            # Add abstract words into list of all words
            words.extend(abstracts.split())
            data_from_count += 1
        except:
            pass

    # calculate the frequency of each word in abstracts
    freq = dict((x, words.count(x)) for x in set(words))

    # = Remove stop words from the list of all words =
    # Read stop words from file
    stop_lines = tuple(open(config.config_dir + "/stopwords", "r"))
    stop_words = []
    for line in stop_lines:
        split = line.split()
        if len(split) > 0 and split[0] != "|" and "|" not in split[0]:
            stop_words.append(split[0])

    # Remove words
    for stp in stop_words:
        freq.pop(stp, None)

    # Remove anything with " in it (Causes javascript problems)
    toRemove = []
    for f in freq:
        if '"' in f:
            toRemove.append(f)

    for f in toRemove:
        freq.pop(f, None)

    i = 0
    print 'Top 5'
    # print a list of sorted frequencies
    with open(config.data_dir + '/abstracts.csv', 'wb') as csvfile:
        abstracts_file = csv.writer(csvfile)
        for w in sorted(freq, key=freq.get, reverse=True):
            if i < 5:
                print w, freq[w]
                i = i+1
            abstracts_file.writerow([w.encode('utf-8'), freq[w]])

    return data_from_count


############################################################
# Try with the authors - these are in a nested dict
############################################################
def authors(papers):

    import hashlib

    authors = []
    for this_paper in papers:
        # Some pmid files dont actually have an authorlist! e.g. 2587412
        # This probably needs to be resolved with pubmed!
        try:
            for this_author in this_paper['author']:
                # There are some entries in the author list that are not actually authors e.g. 21379325 last author
                try:
                    authors.append(this_author['family'])
                    # Create a clean author field. This is the Surname followed by first initial.
                    this_author.update({'clean': this_author['family'] + " " + this_author['given'][0]})
                except:
                    pass

        except:
            logging.warn('No AuthorList for ' + this_paper['IDs']['hash'])

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
            for this_author in this_paper['author']:
                # Create author hash
                hash_object = hashlib.sha256(this_author['clean'])
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
                for con_author in this_paper['author']:

                    if not this_author == con_author:
                        con_hash_object = hashlib.sha256(con_author['clean'])
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
# Output a report on how complete the data is
############################################################
def completeness_report(papers):

    num_papers = len(papers)
    first_authors_count = 0
    for this_paper in papers:
        try:
            if this_paper['author'][0]['family']:
                first_authors_count = first_authors_count + 1
        except:
            next

    print 'Total papers = ' + str(num_papers)
    print 'First authors count = ' + str(first_authors_count)


############################################################
# Try with the FIRST authors - these are in a nested dict
############################################################
def first_authors(papers):

    num_papers = len(papers)
    first_authors = []
    for this_paper in papers:
        try:
            first_authors.append(this_paper['author'][0]['family'])
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
            first_authors_inst.append(this_paper['Extras']['CleanInstitute'])
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

        if 'MedlineCitation' in this_paper:
            if 'MeshHeadingList' in this_paper['MedlineCitation']:
                coverage = coverage + 1

                try:
                    for this_mesh in this_paper['MedlineCitation']['MeshHeadingList']:
                        mesh.append(this_mesh['DescriptorName'])
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


# output a csv file with some info in
def output_csv(papers):

    print '\n###Outputting CSV file###'
    with open(config.data_dir + '/all.csv', 'wb') as csvfile:
        all_file = csv.writer(csvfile)
        for this_paper in papers:
            try:
                title = this_paper['title']
            except:
                title = '???'

            try:
                first_author = this_paper['author'][0]['family']
            except:
                first_author = '???'

            try:
                journal = this_paper['MedlineCitation']['Article']['Journal']['ISOAbbreviation']
            except:
                journal = '???'

            try:
                citations = this_paper['Extras']['Citations']
            except:
                citations = '???'

            try:
                all_file.writerow([title, first_author, journal, citations, this_paper])
            except:
                pass

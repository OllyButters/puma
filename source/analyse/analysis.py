#! /usr/bin/env python

import csv
import logging

############################################################
# Have all the data now, so do something with it


############################################################
# Build a list of all journals and count frequencies of each
def journals(papers):
    print "\n###Journals###"

    num_papers = len(papers)

    journals = []
    for this_paper in papers:
        journals.append(this_paper['MedlineCitation']['Article']['Journal']['ISOAbbreviation'])

    print str(len(journals))+'/'+str(num_papers)
    print str(len(set(journals)))+' different journals'

    # calculate the frequency of each journal
    freq = dict((x, journals.count(x)) for x in set(journals))

    i = 0
    print 'Top 5'

    # print a list of sorted frequencies
    with open('../data/journals.csv', 'wb') as csvfile:
        journals_file = csv.writer(csvfile)
        for w in sorted(freq, key=freq.get, reverse=True):
            if i < 5:
                print w, freq[w]
                i = i+1
            journals_file.writerow([w.encode('utf-8'), freq[w]])


############################################################
# Build a list of Abstracts
def abstracts(pmids, papers):
    print "\n###Abstracts###"

    abstracts = ''
    for this_pmid in pmids:
        try:
            abstracts = abstracts + str(papers[this_pmid]['AbstractText'])
        except:
            pass

    words = abstracts.split()

    # calculate the frequency of each word in abstracts
    freq = dict((x, words.count(x)) for x in set(words))

    i = 0
    print 'Top 5'

    # print a list of sorted frequencies
    with open('../data/abstracts.csv', 'wb') as csvfile:
        abstracts_file = csv.writer(csvfile)
        for w in sorted(freq, key=freq.get, reverse=True):
            if i < 5:
                print w, freq[w]
                i = i+1
            abstracts_file.writerow([w.encode('utf-8'), freq[w]])


############################################################
# Try with the authors - these are in a nested dict
def authors(papers):

    import hashlib

    authors = []
    for this_paper in papers:
        # print papers[this_pmid]
        # Some pmid files dont actually have an authorlist! e.g. 2587412
        # This probably needs to be resolved with pubmed!
        try:
            for this_author in this_paper['author']:
                # There are some entries in the author list that are not actually authors e.g. 21379325 last author
                try:
                    authors.append(this_author['family'])
                except:
                    pass
        except:
            logging.warn('No AuthorList for '+this_paper['IDs']['hash'])

    # print authors
    freq = dict((x, authors.count(x)) for x in set(authors))
    print "\n###Authors###"

    print str(len(set(authors)))+' different authors'
    # print freq

    i = 0
    print 'Top 5'

    with open('../data/authors.csv', 'wb') as csvfile:
        authors_file = csv.writer(csvfile)
        for w in sorted(freq, key=freq.get, reverse=True):
            if i < 5:
                print w, freq[w]
                i = i+1
        # Need to utf-8 encode
            authors_file.writerow([w.encode('utf-8'), freq[w]])

    # === Author Network Object ===
    author_network = {}
    author_network['authors'] = {}
    author_network['connections'] = {}
    for this_paper in papers:
        try:
            for this_author in this_paper['author']:
                # Create author hash
                hash_object = hashlib.sha256(this_author['family'] + this_author['given'])
                author_hash = hash_object.hexdigest()

                # Store author details
                try:
                    author_network['authors'][ author_hash ]
                except:
                    author_network['authors'][ author_hash ] = {}
                    author_network['authors'][ author_hash ]['family'] = this_author['family']
                    author_network['authors'][ author_hash ]['given'] = this_author['given']

                try:
                    author_network['authors'][ author_hash ]['num_papers'] += 1
                except:
                    author_network['authors'][ author_hash ]['num_papers'] = 1

                # Create edges
                for con_author in this_paper['author']:

                    if not this_author == con_author:
                        con_hash_object = hashlib.sha256(con_author['family'] + con_author['given'])
                        con_author_hash = con_hash_object.hexdigest()

                        con_id = ""
                        if author_hash > con_author_hash:
                            con_id = author_hash + "" + con_author_hash
                        else:
                            con_id = con_author_hash + "" + author_hash

                        con_id_hash_object = hashlib.sha256(con_id)
                        con_hash = con_id_hash_object.hexdigest()

                        try:
                            author_network['connections'][ con_hash ]
                        except:
                            author_network['connections'][ con_hash ] = {}
                            author_network['connections'][ con_hash ]['authors'] = []
                            author_network['connections'][ con_hash ]['authors'].append( { 'author_hash': author_hash } )
                            author_network['connections'][ con_hash ]['authors'].append( { 'author_hash': con_author_hash } )

                        try:
                            author_network['connections'][ con_hash ]['num_connections'] += 1
                        except:
                            author_network['connections'][ con_hash ]['num_connections'] = 1

        except:
            pass

    print author_network
    print "\n###Author network###"
    print str(len( author_network['authors'])) + " Authors"
    print str(len( author_network['connections'])) + " Connections"


############################################################
# Try with the FIRST authors - these are in a nested dict
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

    with open('../data/first_authors.csv', 'wb') as csvfile:
        authors_file = csv.writer(csvfile)
        for w in sorted(freq, key=freq.get, reverse=True):
            if i < 5:
                print w, freq[w]
                i = i+1
                # Need to utf-8 encode
            authors_file.writerow([w.encode('utf-8'), freq[w]])


############################################################
# Try with the FIRST authors INSTITUTE- these are in a nested dict
def inst(papers):
    num_papers = len(papers)

    first_authors_inst = []
    for this_paper in papers:
        try:
            first_authors_inst.append(this_paper['Extras']['CleanInstitute'])
        except:
            pass

    # print authors
    freq = dict((x, first_authors_inst.count(x)) for x in set(first_authors_inst))
    print "\n###First authors institute###"

    print str(len(first_authors_inst))+'/'+str(num_papers)
    print str(len(set(first_authors_inst)))+' different first author institutes'
    # print str(not_matched)+' institutes not matched with lookup'
    # print freq

    i = 0
    print 'Top 5'

    with open('../data/first_authors_inst.csv', 'wb') as csvfile:
        authors_file = csv.writer(csvfile)
        for w in sorted(freq, key=freq.get, reverse=True):
            if i < 5:
                print w, freq[w]
                i = i+1
        # Need to utf-8 encode
            authors_file.writerow([w.encode('utf-8'), freq[w]])


############################################################
# Try with the mesh headings - these are in a nested dict
def mesh(papers):
    num_papers = len(papers)

    mesh = []
    coverage = 0
    for this_paper in papers:

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

    with open('../data/mesh.csv', 'wb') as csvfile:
        mesh_file = csv.writer(csvfile)
        for w in sorted(freq, key=freq.get, reverse=True):
            if i < 5:
                print w, freq[w]
                i = i+1
            # Need to utf-8 encode
            mesh_file.writerow([w.encode('utf-8'), freq[w]])


# output a csv file with some info in
def output_csv(papers):

    print '\n###Outputting CSV file###\n'
    with open('../data/all.csv', 'wb') as csvfile:
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

            # try:
            #    year         = papers[0]['Year']
            # except:
            #    year = '???'

            # try:
            #    month        = papers[0]['Month']
            # except:
            #    month = '???'

            try:
                citations = this_paper['Extras']['Citations']
            except:
                citations = '???'

            try:
                all_file.writerow([title, first_author, journal, citations, this_paper])
            except:
                print 'Failing on '+str(this_paper['IDs']['hash'])

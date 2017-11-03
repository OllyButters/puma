#! /usr/bin/env python

import csv
import logging
import config.config as config
import os
import shutil

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
# Build a list of words used in Abstracts
############################################################
# def abstracts(papers):
def DEPRICATED(papers):
    print "\n###Abstracts###"

    words = []
    data_from_count = 0
    # Go through all papers
    for this_paper in papers:
        try:
            # Get abstract text
            abstracts = str(this_paper['clean']['abstract'])

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
# Build a list of words used in 'item'
############################################################
def word_frequencies(papers, item):
    print("\n###" + item + "###")

    all_words = []
    data_from_count = 0
    # Go through all papers
    for this_paper in papers:
        try:
            # Get abstract text
            text = str(this_paper['clean']['item'])

            # Remove punctuation and esacpe characters that will cause a problem
            text = text.lower()
            text = text.replace(",", " ")
            text = text.replace(".", " ")
            text = text.replace(":", " ")
            text = text.replace(";", " ")
            text = text.replace("'", "\'")
            text = text.replace('"', ' ')

            # Add abstract words into list of all words
            all_words.extend(text.split())
            data_from_count += 1
        except:
            pass

    # calculate the frequency of each word in abstracts
    freq = dict((x, all_words.count(x)) for x in set(all_words))

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
        freq.pop(stp, None)

    i = 0
    print 'Top 5'

    # Output the data to file and print the top 5 to screen.
    with open(config.data_dir + '/' + item + '.csv', 'wb') as csvfile:
        output_file = csv.writer(csvfile)
        for w in sorted(freq, key=freq.get, reverse=True):
            if i < 5:
                print w, freq[w]
                i = i+1
            output_file.writerow([w.encode('utf-8'), freq[w]])

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
# Build an HTML report of the status of the importat fields
################################################################################
def coverage_report(papers):

    print 'Building coverage report'

    cov_css = '''
                .missing_required {background-color: red;}
                .missing_good_to_have {background-color: orange;}
                tr:nth-child(even) {background-color: #f2f2f2}
                th {background-color: #4CAF50; color: white;}
                td, th {padding: 0.2em;}
                a:visited {color: red;}
    '''

    cov_html = '<table>'
    cov_html += '''<tr>
                    <th>Hash</th>
                    <th>Zotero</th>
                    <th>DOI</th>
                    <th>PMID</th>
                    <th>PMID DOI<br/>Lookup</th>
                    <th>PMID title<br/>Lookup</th>
                    <th>Title</th>
                    <th>Keywords</th>
                    <th>Abstract</th>
                    <th>First<br/>Author</th>
                    <th>First<br/>Author<br/>affil</th>
                    <th>Clean<br/>Inst</th>
                    <th>Clean<br/>Date</th>
                    <th>Journal</th>
                    <th>Scopus<br/>Citations</th>
                </tr>'''

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
    status['clean_date'] = 0
    status['journal'] = 0
    status['scopus'] = 0

    for this_paper in papers:
        cov_html += '<tr>'

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
            cov_html += '<td class="missing_good_to_have">???</td>'

        # Clean date
        try:
            clean_date = this_paper['clean']['clean_date']['year']
            if clean_date != '':
                cov_html += '<td>OK</td>'
                status['clean_date'] = status['clean_date'] + 1
            else:
                raise Exception()
        except:
            cov_html += '<td class="missing_required">???</td>'

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
            scopus = this_paper['clean']['citations']['scopus']['count']
            if scopus != '':
                cov_html += '<td>OK</td>'
                status['scopus'] = status['scopus'] + 1
            else:
                raise Exception()
        except:
            cov_html += '<td class="missing_good_to_have">???</td>'

        cov_html += '</tr>'
    cov_html += '</table>'

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
    status_table += '<td>' + str(int(round(100*status['clean_date']/number_of_papers))) + '</td>'
    status_table += '<td>' + str(int(round(100*status['journal']/number_of_papers))) + '</td>'
    status_table += '<td>' + str(int(round(100*status['scopus']/number_of_papers))) + '</td>'
    status_table += '</tr></table>'

    # Title
    title = '<h1>' + config.project_details['short_name'] + '</h1>'

    output_text = '<html><head><style>' + cov_css + '</style></head><body>' + title + status_table + '<br/><br/>' + cov_html + '</body></html>'
    coverage_file = open(config.html_dir + '/coverage_report.html', 'w')
    print >> coverage_file, output_text

    # put a copy of all the processed files in the web tree
    if os.path.exists(config.html_dir + '/status/cleaned'):
        shutil.rmtree(config.html_dir + '/status/cleaned')
    shutil.copytree(config.cache_dir + '/processed/cleaned', config.html_dir + '/status/cleaned')

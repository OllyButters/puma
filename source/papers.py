#! /usr/bin/env python

#12/9/14

#Seems that there are a few that don't have authors, so the lastname call fails. eg 24770850
#Current batch has a non ASCII character.

#todo
#wordcloud of abstracts
#Average position in list of authors
#institutions - start to put rules in place - will take a while!
#will need to deal with missing data somehow.
#seems to be lots of missing mesh headings
#Google docs has an API, so could pull all the pmids from that. Would be a simple interface.
#https://developers.google.com/google-apps/spreadsheets/#working_with_list-based_feeds

import csv
import json
import os.path

import tools

from Bio import Entrez
Entrez.email = "olly.butters@bristol.ac.uk"     # Always tell NCBI who you are


###########################################################
#list of PMIDs - this would be fed from somewhere
#three of pauls
#pmids=['25085103','15831561','16154023']

#last ten alspac
#pmids=['24963150','24924479','24930394','24952709','24945404','23841856','24848214','23895510','24158349']
#'24927274' This one doesnt work

#Read in the list of PMIDs from an external csv file in this directory.
#The file must be just pmids - one per line, and no extra lines at the end.
pmids = []
with open('../pmids.csv', 'rb') as csvfile:
    f = csv.reader(csvfile)
    for row in f:
        #print row[0]
        pmids.append(row[0])

#print pmids

num_pmids=len(pmids)
print str(num_pmids)+' PMIDs to process'


###########################################################
#Try to build a papers dictionary with PMID as the index
papers = {}

for this_pmid in pmids:

    print 'Working on '+this_pmid


    #Using efetch

    #Build a cache of all the pmid data so we don't keep downloading it.
    #Check that cache first when looking for a PMID, if it's not there then
    #go and get the data.
    file_name='../cache/'+this_pmid 
    if not os.path.isfile(file_name):
        #This PMID data not in cache, so download it.
        handle = Entrez.efetch(db="pubmed", id=this_pmid, retmode="xml")
        record = Entrez.read(handle)
        
        #Should check to see if anything sensible was returned - eg check size

        #Save it for later
        file_name='../cache/'+this_pmid
        fo = open(file_name, 'wb')
        fo.write(json.dumps(record[0], indent=4))
        fo.close()

    #Open and parse cached file
    fo = open(file_name, 'r')
    record = json.loads(fo.read())
    fo.close()


    #print all
    #print '#####ALL'
    #print record


    #print '#####Authors'
    #print all authors
    #print record[0]['MedlineCitation']['Article']['AuthorList']
    
    
    #print out lastname and affiliation for each author
    #for this_one in record[0]['MedlineCitation']['Article']['AuthorList']:
     #   print this_one['LastName']
      #  print this_one['Affiliation']
               
    #all the info for this paper
    this_paper = {}
    this_paper['pmid']=this_pmid
    this_paper['ArticleTitle'] = record['MedlineCitation']['Article']['ArticleTitle']
    this_paper['Journal'] = record['MedlineCitation']['Article']['Journal']['ISOAbbreviation']
    
    #Mesh keywords
    try:
        this_paper['MeshHeadingList'] = record['MedlineCitation']['MeshHeadingList']
    except:
        #Should log this
        print 'No mesh headings'

    #Author list
    try:
        this_paper['AuthorList'] = record['MedlineCitation']['Article']['AuthorList']
    except:
        print 'No Authors listed!'

    
    #Add this_paper info to the main dict
    papers[this_pmid]=this_paper


#Save it for later
file_name='../data/summary'
fo = open(file_name, 'wb')
fo.write(json.dumps(papers, indent=4))
fo.close()

print '####All papers'
#print papers

print str(num_pmids)+' PMIDs to process'
print str(len(papers))+' papers processed'


############################################################
#Have all the data now, so do something with it

############################################################
#Build a list of all journals and count frequencies of each
print "\n###Journals###"

journals = []
for this_pmid in pmids:
    journals.append(papers[this_pmid]['Journal'])

print str(len(journals))+'/'+str(num_pmids)
print str(len(set(journals)))+' different journals'

#calculate the frequency of each journal
freq = dict((x,journals.count(x)) for x in set(journals))

i=0
print 'Top 5'

#print a list of sorted frequencies
with open('../data/journals.csv', 'wb') as csvfile:
    journals_file = csv.writer(csvfile)
    for w in sorted(freq, key=freq.get, reverse=True):
        if i<5:
            print w, freq[w]
            i=i+1
        journals_file.writerow([w.encode('utf-8'), freq[w]])





############################################################
#Try with the authors - these are in a nested dict
authors = []
for this_pmid in pmids:
    #print this_pmid
    for this_author in papers[this_pmid]['AuthorList']:
        #There are some entries in the author list that are not actually authors e.g. 21379325 last author
        try:
            authors.append(this_author['LastName'])
        except:
            pass


#print authors
freq = dict((x,authors.count(x)) for x in set(authors))
print "\n###Authors###"

print str(len(set(authors)))+' different authors'
#print freq

i=0
print 'Top 5'

with open('../data/authors.csv', 'wb') as csvfile:
    authors_file = csv.writer(csvfile)
    for w in sorted(freq, key=freq.get, reverse=True):
        if i<5:
            print w, freq[w]
            i=i+1
        #Need to utf-8 encode
        authors_file.writerow([w.encode('utf-8'), freq[w]])


############################################################
#Try with the FIRST authors - these are in a nested dict
first_authors = []
for this_pmid in pmids:
    first_authors.append(papers[this_pmid]['AuthorList'][0]['LastName'])

#print authors
freq = dict((x,first_authors.count(x)) for x in set(first_authors))
print "\n###First authors###"

print str(len(first_authors))+'/'+str(num_pmids)
print str(len(set(first_authors)))+' different first authors'

#print freq

i=0
print 'Top 5'

with open('../data/first_authors.csv', 'wb') as csvfile:
    authors_file = csv.writer(csvfile)
    for w in sorted(freq, key=freq.get, reverse=True):
        if i<5:
            print w, freq[w]
            i=i+1
        #Need to utf-8 encode
        authors_file.writerow([w.encode('utf-8'), freq[w]])

############################################################
#Try with the FIRST authors INSTITUTE- these are in a nested dict
first_authors_inst = []
for this_pmid in pmids:
    try:
        first_authors_inst.append(papers[this_pmid]['AuthorList'][0]['Affiliation'])
    except:
        pass

tools.clean_institution(first_authors_inst)

#print authors
freq = dict((x,first_authors_inst.count(x)) for x in set(first_authors_inst))
print "\n###First authors institute###"

print str(len(first_authors_inst))+'/'+str(num_pmids)
print str(len(set(first_authors_inst)))+' different first author institutes'
#print freq

i=0
print 'Top 5'

with open('../data/first_authors_inst.csv', 'wb') as csvfile:
    authors_file = csv.writer(csvfile)
    for w in sorted(freq, key=freq.get, reverse=True):
        if i<5:
            print w, freq[w]
            i=i+1
        #Need to utf-8 encode
        authors_file.writerow([w.encode('utf-8'), freq[w]])



############################################################
#Try with the mesh headings - these are in a nested dict
mesh = []
coverage=0
for this_pmid in pmids:
    #print this_pmid
    
    if 'MeshHeadingList' in papers[this_pmid]:
        coverage = coverage + 1
    
    try:
        for this_mesh in papers[this_pmid]['MeshHeadingList']:
            mesh.append(this_mesh['DescriptorName'])
    except:
        #Do nothing
        pass

#print mesh
freq = dict((x,mesh.count(x)) for x in set(mesh))
print "\n###Mesh###"

#print freq

print str(coverage)+'/'+str(num_pmids)
print str(len(set(mesh)))+' different mesh headings'

i=0
print 'Top 5'

with open('../data/mesh.csv', 'wb') as csvfile:
    mesh_file = csv.writer(csvfile)
    for w in sorted(freq, key=freq.get, reverse=True):
        if i<5:
            print w, freq[w]
            i=i+1
        #Need to utf-8 encode
        mesh_file.writerow([w.encode('utf-8'), freq[w]])

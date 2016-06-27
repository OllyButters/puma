#! /usr/bin/env python

#Have a go at tidying up the mess that is first author institution.
#Essentially go through each institution and see if it matches a patten
#in the institute_cleaning.csv file. If it does then replace it with a
#standard name.
def clean_institution(paper_list):
    import csv
    import re
    import logging
    import json

    logging.info('Starting institute cleaning')

    #Read in config file
    pattern = []
    replacements = []
    with open('../config/institute_cleaning.csv', 'rb') as csvfile:
        f = csv.reader(csvfile)
        for row in f:
            try:
                #Check it is not a comment string first.
                if(re.match('#',row[0])):
                    continue

                #Check for blank lines
                if row[0] == '':
                    continue

                #If there is a second element in this row then carry on
                pattern.append(row[0])
                replacements.append(row[1])
            except:
                pass

    #print pattern
    #print replacements

    logging.info('Config read in, starting processing')

    #Cycle through institute checking the whole substitution list.
    #Stop when the first one matches.
    number_not_matched=0
    for this_paper in paper_list:

        #open the file and parse it
        file_name='../cache/raw/'+this_paper
        print file_name

        with open(file_name) as fo:
            papers=json.load(fo)

        try:
            #institute = papers[this_pmid]['AuthorList'][0]['Affiliation']
            #institute = papers[this_pmid]['AuthorList'][0]['AffiliationInfo'][0]['Affiliation']
            #institute = papers[this_pmid]['AuthorList'][0]['AffiliationInfo'][0]['Affiliation']
            #print papers
            print '============='
            print papers[0]['author'][0]['affiliation'][0]['name']
            institute = papers[0]['author'][0]['affiliation'][0]['name']
            #exit()

        except:
            logging.warn('Could not find an affiliation for %s', this_paper)
            continue

        for y in range(0,len(pattern)):
            logging.debug('%s %s %s', institute, pattern[y], replacements[y])
            temp = re.search(pattern[y], institute, re.IGNORECASE)
            if(temp>0):
                logging.info('ID:%s. %s MATCHES %s REPLACEDBY %s', this_paper, institute, pattern[y], replacements[y])
                #papers[this_paper]['Extras']['CleanInstitute'] = replacements[y]
                papers[0]['Extras'] = {}
                papers[0]['Extras']['CleanInstitute'] = replacements[y]
                break

            if(y==len(pattern)-1):
                logging.info('No match for %s. ID:%s', institute, this_paper)
                logging.warn('No match for %s. ID:%s', institute, this_paper)
                number_not_matched+=1

        #Save it for later
        file_name='../cache/processed/'+this_paper
        fo = open(file_name, 'wb')
        fo.write(json.dumps(papers, indent=4))
        fo.close()


    print 'Cleaning institutions'
    print str(len(paper_list)-number_not_matched)+'/'+str(len(paper_list))+' cleaned'


    return number_not_matched



#Go through the deltas directory and apply any changes that are needed
def do_deltas(papers):
    import csv
    import json
    import logging
    import os

    delta_dir = '../config/deltas/'


    deltas = os.listdir(delta_dir)

    print deltas

    for this_delta in deltas:
        #delta_file = '8680184'

        delta_path = delta_dir+this_delta

        fo = open(delta_path, 'r')
        record = json.loads(fo.read())
        fo.close()

        try:
            papers[this_delta]['Year'] = record['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['Year']
        except:
            print 'FAIL'

#! /usr/bin/env python

#Have a go at tidying up the mess that is first author institution.
#Essentially go through each institution and see if it matches a patten
#in the institut_cleaning.csv file. If it does then replace it with a 
#standard name.
def clean_institution(pmids,papers):
    import csv
    import re
    import logging

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

                #If there is a second element in this row then carry on
                pattern.append(row[0])
                replacements.append(row[1])
            except:
                pass

    #Cycle through institute checking the whole substitution list.
    #Stop when the first one matches.
    number_not_matched=0
    for this_pmid in pmids:
        
        try:
            institute = papers[this_pmid]['AuthorList'][0]['Affiliation']
        except:
            continue
        
        for y in range(0,len(pattern)):
            #print pattern[y]
            temp = re.search(pattern[y], institute)
            if(temp>0):
                logging.info('%s MATCHES %s REPLACEDBY %s', institute, pattern[y], replacements[y])
                papers[this_pmid]['Extras']['CleanInstitute'] = replacements[y]
                break
            
            if(y==len(pattern)-1):
                logging.info('No match for %s', institute)
                logging.warn('No match for %s', institute)
                number_not_matched+=1

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

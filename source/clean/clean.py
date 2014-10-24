#! /usr/bin/env python

#Have a go at tidying up the mess that is first author institution.
#Essentially go through each institution and see if it matches a patten
#in the institut_cleaning.csv file. If it does then replace it with a 
#standard name.
#def clean_institution_course2(inst):
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
                #print inst[x].encode('latin-1')+' ##matches## '+pattern[y]+' ##replaced by## '+replacements[y]
                logging.info('%s MATCHES %s REPLACEDBY %s', institute, pattern[y], replacements[y])
                papers[this_pmid]['Extras']['CleanInstitute'] = replacements[y]
                break
            
            if(y==len(pattern)-1):
                #print 'No match for '+inst[x].encode('latin-1')           
                logging.info('No match for %s', institute)
                logging.warn('No match for %s', institute)
                number_not_matched+=1

    return number_not_matched
#    print str(number_not_matched)+' not matched with one in institution lookup file'



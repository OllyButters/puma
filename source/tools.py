#! /usr/bin/env python

#Have a go at tidying up the mess that is first author institution.
#Essentially go through each institution and see if it matches a patten
#in the institut_cleaning.csv file. If it does then replace it with a 
#standard name.
def clean_institution_course2(inst):
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
    for x in range(0, len(inst)):
        for y in range(0,len(pattern)):
            #print pattern[y]
            temp = re.search(pattern[y], inst[x])
            if(temp>0):
                #print inst[x].encode('latin-1')+' ##matches## '+pattern[y]+' ##replaced by## '+replacements[y]
                logging.info('%s MATCHES %s REPLACEDBY %s', inst[x], pattern[y], replacements[y])
                inst[x] = replacements[y]
                break
            
            if(y==len(pattern)-1):
                #print 'No match for '+inst[x].encode('latin-1')           
                logging.info('No match for %s', inst[x])
                logging.warn('No match for %s', inst[x])
                number_not_matched+=1

    return number_not_matched
#    print str(number_not_matched)+' not matched with one in institution lookup file'


#This approach tries to clean the junk out of individual institutions. I fear there is
#too much variation.
def clean_institution_fine(inst):
    import re
    for x in range(0, len(inst)):
        #print inst[x]
        
        #delete email address
        inst[x] = re.sub(r'[\w\.-]+@[\w\.-]+', r'', inst[x])

        #delete 'Electronic address:'
        inst[x] = re.sub(r'Electronic address:', r'', inst[x])

        #replace United Kingdom with UK
        inst[x] = re.sub(r'United Kingdom', r'UK', inst[x])

        #replace & with and
        inst[x] = re.sub(r'&', r'and', inst[x])

        #get rid of Oakfield house
        inst[x] = re.sub(r'BS8 2BN,? ?', r'', inst[x])
        inst[x] = re.sub(r'BS82BN,? ?', r'', inst[x])
        inst[x] = re.sub(r'Oakfield House, Oakfield Grove, Clifton, ', r'', inst[x])
        inst[x] = re.sub(r'Oakfield House, Oakfield Grove, ', r'', inst[x])

        inst[x] = re.sub(r'Canynge Hall, 39 Whatley Road, ', r'', inst[x])

        inst[x] = re.sub(r'BS8 2PS, ', r'', inst[x])

        inst[x] = re.sub(r'University of Bristol, Bristol,? ?', r'University of Bristol, ', inst[x])
        #inst[x] = re.sub(r'University of Bristol, Bristol ', r'University of Bristol, ', inst[x])

        inst[x] = re.sub(r'Medical Research Council', r'MRC', inst[x])

        inst[x] = re.sub(r' at the', r',', inst[x])

        inst[x] = re.sub(r'Avon,', r'', inst[x])


        #Get rid of any trailing whitespace
        inst[x] = inst[x].rstrip()

        #print inst[x]

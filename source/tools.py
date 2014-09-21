#! /usr/bin/env python

#Have a go at tidying up the mess that is first author institution.

def clean_institution_course(inst):
    import re

    for x in range(0, len(inst)):
        #print inst[x]
        
        #Do email addresses first - my guess is that people are more likely to 
        #have only one email address than multiple institution addresses...
        inst[x] = re.sub(r'bbk.ac.uk', r'Birbeck', inst[x])
        inst[x] = re.sub(r'(?i)bristol.ac.uk|bris.ac.uk', r'University of Bristol', inst[x])
        inst[x] = re.sub(r'bpmf.ac.uk', r'bpmf.ac.uk\?', inst[x])
        inst[x] = re.sub(r'cam.ac.uk', r'University of Cambridge', inst[x])
        inst[x] = re.sub(r'city.ac.uk', r'City University', inst[x])
        inst[x] = re.sub(r'ex.ac.uk', r'University of Essex', inst[x])
        inst[x] = re.sub(r'einstein.yu.edu', r'Albert Einstein College of Medicine', inst[x])
        inst[x] = re.sub(r'gla.ac.uk', r'University of Glasgow', inst[x])
        inst[x] = re.sub(r'harvard.edu', r'Harvard University', inst[x])
        inst[x] = re.sub(r'jhsph.edu|jhmi.edu', r'Johns Hopkins University', inst[x])
        inst[x] = re.sub(r'kcl.ac.uk', r'KCL', inst[x])
        inst[x] = re.sub(r'liverpool.ac.uk', r'University of Liverpool', inst[x])
        inst[x] = re.sub(r'(?i)lshtm.ac.uk', r'LSHTM', inst[x])
        inst[x] = re.sub(r'nottingham.ac.uk', r'University of Nottingham', inst[x])
        inst[x] = re.sub(r'ox.ac.uk', r'University of Oxford', inst[x])
        inst[x] = re.sub(r'pds.ac.uk', r'Peninsula College of Medicine & Dentistry', inst[x])
        inst[x] = re.sub(r'qmw.ac.uk|qmul', r'Queen Mary University London', inst[x])
        inst[x] = re.sub(r'qub.ac.uk', r'Queens University', inst[x])
        inst[x] = re.sub(r'sanger.ac.uk', r'Sanger Institute', inst[x])
        inst[x] = re.sub(r'soton.ac.uk', r'University of Southampton', inst[x])
        inst[x] = re.sub(r'ucl.ac.uk', r'UCL', inst[x])
        inst[x] = re.sub(r'uwe.ac.uk', r'UWE', inst[x])


        #Do the UK institution addresses.
        inst[x] = re.sub(r'^.*University of Aberdeen.*$', r'University of Aberdeen', inst[x])
        inst[x] = re.sub(r'^.*University of Bath.*$', r'University of Bath', inst[x])
        inst[x] = re.sub(r'^.Birbeckn.*$', r'Birbeck', inst[x])
        inst[x] = re.sub(r'^.*University of Birmingham.*$|^.*Birmingham University.*$', r'University of Birmingham', inst[x])
        inst[x] = re.sub(r'^.*University of Brighton.*$', r'University of Brighton', inst[x])
        inst[x] = re.sub(r'^.*University of Bristol.*$', r'University of Bristol', inst[x])
        inst[x] = re.sub(r'^.*Brunel University.*$', r'Brunel University', inst[x])
        inst[x] = re.sub(r'^.*University of Cambridge.*$', r'University of Cambridge', inst[x])
        inst[x] = re.sub(r'^.*Cardiff University.*$', r'Cardiff University', inst[x])
        inst[x] = re.sub(r'^.*City University.*$', r'City University', inst[x])
        inst[x] = re.sub(r'^.*University of Coventry.*$|^.*Coventry University.*$', r'University of Coventry', inst[x])
        inst[x] = re.sub(r'^.*DeMontfort University.*$', r'DeMontfort University', inst[x])
        inst[x] = re.sub(r'^.*University of Dundee.*$', r'University of Dundee', inst[x])
        inst[x] = re.sub(r'^.*University of Durham.*$', r'University of Durham', inst[x])
        inst[x] = re.sub(r'^.*University of Edinburgh.*$', r'University of Edinburgh', inst[x])
        inst[x] = re.sub(r'^.*University of Essex.*$', r'University of Essex', inst[x])
        inst[x] = re.sub(r'^.*University of Exeter.*$', r'University of Exeter', inst[x])
        inst[x] = re.sub(r'^.*University of Glasgow.*$', r'University of Glasgow', inst[x])
        inst[x] = re.sub(r'^.*Imperial.*$', r'Imperial', inst[x])
        inst[x] = re.sub(r'^.*Institute of Education.*$', r'IoE', inst[x])
        inst[x] = re.sub(r'^.*King\'?s College,? London.*$', r'KCL', inst[x])
        inst[x] = re.sub(r'^.*Lancaster University.*$', r'Lancaster University', inst[x])
        inst[x] = re.sub(r'^.*University of Leeds.*$', r'University of Leeds', inst[x])
        inst[x] = re.sub(r'^.*University of Leicester.*$', r'University of Leicester', inst[x])
        inst[x] = re.sub(r'^.*University of Liverpool.*$', r'University of Liverpool', inst[x])
        inst[x] = re.sub(r'^.*London Metropolitan University.*$', r'London Metropolitan University', inst[x])
        inst[x] = re.sub(r'^.*University of Manchester.*$', r'University of Manchester', inst[x])
        inst[x] = re.sub(r'^.*University of Munich.*$', r'Universtiy of Munich', inst[x])
        inst[x] = re.sub(r'^.*Newcastle University.*$', r'Newcastle University', inst[x])
        inst[x] = re.sub(r'^.*University of Newcastle.*$', r'Newcastle University', inst[x])
        inst[x] = re.sub(r'^.*University of Nottingham.*$', r'University of Nottingham', inst[x])
        inst[x] = re.sub(r'^.*University of Oxford.*$', r'University of Oxford', inst[x])
        inst[x] = re.sub(r'^.*Peninsula.*$', r'Peninsula College of Medicine & Dentistry', inst[x])
        inst[x] = re.sub(r'^.*Plymouth.*$', r'Plymouth University', inst[x])
        inst[x] = re.sub(r'^.*Queen\'s University.*$', r'Queens University', inst[x])
        inst[x] = re.sub(r'^.*University of Sheffield.*$', r'University of Sheffield', inst[x])
        inst[x] = re.sub(r'^.*University of Southampton.*$', r'University of Southampton', inst[x])
        inst[x] = re.sub(r'^.*University of Strathclyde.*$', r'University of Strathclyde', inst[x])
        inst[x] = re.sub(r'^.*University of Stirling.*$', r'University of Stirling', inst[x])
        inst[x] = re.sub(r'^.*University of Surrey.*$', r'University of Surrey', inst[x])
        inst[x] = re.sub(r'^.*UCL.*$|^.*University College London.*$', r'UCL', inst[x])
        inst[x] = re.sub(r'^.*University of East Anglia.*$', r'UEA', inst[x])
        inst[x] = re.sub(r'^.*University of Warwick.*$', r'University of Warwick', inst[x])
        inst[x] = re.sub(r'^.*University of the West of England.*$|uwe.ac.uk', r'UWE', inst[x])
        inst[x] = re.sub(r'^.*University of York.*$', r'University of York', inst[x])

        #International
        inst[x] = re.sub(r'^.*23andMe.*$', r'23andMe', inst[x])
        inst[x] = re.sub(r'^.*University of Adelaide.*$', r'University of Adelaide', inst[x])
        inst[x] = re.sub(r'^.*University of Alabama.*$', r'University of Alabama', inst[x])
        inst[x] = re.sub(r'^.*University of Amsterdam.*$', r'University of Amsterdam', inst[x])
        inst[x] = re.sub(r'^.*Boston University.*$', r'Boston University', inst[x])
        inst[x] = re.sub(r'^.*University of California.*$', r'University of California', inst[x])
        inst[x] = re.sub(r'^.*University of Chicago.*$', r'University of Chicago', inst[x])
        inst[x] = re.sub(r'^.*Emory University.*$', r'Emory University', inst[x])
        inst[x] = re.sub(r'^.*University of Groningen.*$', r'University of Gronigen', inst[x])
        inst[x] = re.sub(r'^.*University of Hong Kong.*$', r'University of Hong Kong', inst[x])
        inst[x] = re.sub(r'^.*Hong Kong Polytechnic University.*$', r'Hong Kong Polytechnic University', inst[x])
        inst[x] = re.sub(r'^.*Harokopio University.*$', r'Harokopio University', inst[x])
        inst[x] = re.sub(r'^.*Harvard University.*$', r'Harvard University', inst[x])
        inst[x] = re.sub(r'^.*Indiana University.*$', r'Indiana University', inst[x])
        inst[x] = re.sub(r'^.*University of Iowa.*$', r'University of Iowa', inst[x])
        inst[x] = re.sub(r'^.*McGill University.*$', r'McGill University', inst[x])
        inst[x] = re.sub(r'^.*University of Melbourne.*$', r'University of Melbourne', inst[x])
        inst[x] = re.sub(r'^.*University of North Carolina.*$', r'University of North Carolina', inst[x])
        inst[x] = re.sub(r'^.*NYU.*$', r'NYU', inst[x])
        inst[x] = re.sub(r'^.*University of Oregon.*$', r'University of Oregon', inst[x])
        inst[x] = re.sub(r'^.*University of Queensland.*$', r'University of Queensland', inst[x])
        inst[x] = re.sub(r'^.*University of Rochester.*$', r'University of Rochester', inst[x])
        inst[x] = re.sub(r'^.*University of Saskatchewan.*$', r'University of Saskatchewan', inst[x])
        inst[x] = re.sub(r'^.*University of South Australia.*$', r'University of South Australia', inst[x])
        inst[x] = re.sub(r'^.*University of South Carolina.*$', r'University of South Carolina', inst[x])
        inst[x] = re.sub(r'^.*University of Sydney.*$', r'University of Sydney', inst[x])
        inst[x] = re.sub(r'^.*Syracuse University.*$', r'Syracuse University', inst[x])
        inst[x] = re.sub(r'^.*University of Toronto.*$', r'University of Toronto', inst[x])
        inst[x] = re.sub(r'^.*Virginia Commonwealth University.*$', r'Virginia Commonwealth University', inst[x])
        inst[x] = re.sub(r'^.*Yale University.*$', r'Yale University', inst[x])
        inst[x] = re.sub(r'^.*VU University.*$', r'VU University', inst[x])

        #Some special ALSPAC munging
        inst[x] = re.sub(r'^.*Oakfield Grove.*$', r'University of Bristol', inst[x])
        inst[x] = re.sub(r'^.*Canynge Hall.*$', r'University of Bristol', inst[x])

        #Anything left that has bristol in it
        #inst[x] = re.sub(r'^.*Bristol.*$', r'Bristol MISC', inst[x])


#This approach tries to clean individual institutions. I fear there is
#too much variation.
def clean_institution_fine(inst):
    import re
    #print "hello"
    #print inst[0]
    #    for this_inst in inst:
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

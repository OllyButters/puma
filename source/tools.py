#! /usr/bin/env python

def clean_institution(inst):
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

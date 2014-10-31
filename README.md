#Overview#

The general point of this 'pipeline' is to get meta-data about published papers and do some analysis/plotting of it. The general worflow is:

* Start with a list of papers (e.g. with a PubMed ID), and get the meta-data for each.
* Pull out the required meta-data from the originals into a cut down version.
* Clean it.
* Analyse it.
* Make some web pages showing it.
* Build a bibliographic reference of all of them for redistribution.


#Pre-requisites#

##Dependencies##

###Biopython###
To pull out any PubMed data the biopython library is needed:

sudo apt-get install python-biopython

main docs here
http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec%3Aefetch

http://www.nlm.nih.gov/bsd/licensee/elements_descriptions.html

maybe useful info on field descriptions
http://www.nlm.nih.gov/bsd/mms/medlineelements.html

###Google docs###
Can interact with google docs - not currently using this.

sudo apt-get install python-gdata

##Source data##

Currently the whole process kicks off with a list of PubMed IDs. This is generated from a google docs spreadsheet and downloaded to the inputs directory.


#Program flow#

Run the program from /source/papers.py

##Get##


##Clean##

##Add##

##Analyse##

##HTML##


#Notes#
##Dates##
How do you define when a paper was published? Is it the day it is accepted for publication, the day an e-print is published online, or the day the physical paper copy is released? PubMed have thought about this and defined some rules - 
http://www.nlm.nih.gov/bsd/licensee/elements_article_source.html
Essentially they are saying it is up to the publisher to decide when the date should be. An example is e.g. 
http://www.ncbi.nlm.nih.gov/pubmed/20860432
The e-pub was in 2010, but the paper version was in 2011. The publishers (and hence PubMed) go for the later date in this case. 


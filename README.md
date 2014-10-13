
to install dependencies:


=Biopython
sudo apt-get install python-biopython

main docs here
http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec%3Aefetch

=Google docs
sudo apt-get install python-gdata



maybe useful info on field descriptions
http://www.nlm.nih.gov/bsd/mms/medlineelements.html

The general workflow is:
- Fill in the google doc spreadsheet
- Export the PMIDs from there to a csv file
- Harvest the Pubmed data, making a cache of it as we go
- Read in the data from the cache, munging it as we go
- Do the analysis on the munged data

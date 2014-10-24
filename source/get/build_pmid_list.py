#! /usr/bin/env python

###############################################
#Get the pmids from google docs and parse them
#into a single list
###############################################

#19/9/14

import urllib
import csv

#Get the URL of shared spreadsheet and put it here. Note the /export?format=csv bit on the end
url = 'https://docs.google.com/spreadsheets/d/1w2yUIbzSqlnswFA40kcB-zV3ORfwkQpP7Nd53l8Zuhc/export?format=csv'
filename = '../inputs/google_docs.csv'

urllib.urlretrieve(url, filename)
urllib.urlretrieve(url, filename)


#Now we need to do somthing with it!
pmids = []
with open('../inputs/google_docs.csv', 'rb') as csvfile:
    f = csv.reader(csvfile)
    for row in f:
        print row[0]
        #pmids.append(row[0])


#This sheet needs to be published to the web
from gdata.spreadsheet.service import SpreadsheetsService
#key = '0Aip8Kl9b7wdidFBzRGpEZkhoUlVPaEg2X0F2YWtwYkE'
key = '1w2yUIbzSqlnswFA40kcB-zV3ORfwkQpP7Nd53l8Zuhc'


client = SpreadsheetsService()
#feed = client.GetWorksheetsFeed(key, visibility='public', projection='basic')

feed = client.GetCellsFeed(key, visibility='public', projection='full')

for sheet in feed.entry:
  print sheet.title.text
  print sheet.content.text
  #print sheet


#feed = client.GetMedia()



#from apiclient import errors
# ...

#download_url = drive_file.get(key)
#if download_url:
#    resp, content = service._http.request(download_url)
#    if resp.status == 200:
#        print 'Status: %s' % resp
#        #return content
#    else:
#        print 'An error occurred: %s' % resp
#        #return None
#else:
#    # The file doesn't have any content stored on Drive.
#    #return None
#    print "hi"

#import gdata.docs.service


#gdata.docs.service.DownloadResource()

#exit
# Create a client class which will make HTTP requests with Google Docs server.
#client = gdata.docs.service.DocsService()


#worksheet = spreadsheet.sheet1
#contents = []
#for rows in worksheet.get_all_values():
#    contents.append(rows)

# Authenticate using your Google Docs email address and password.
#client.ClientLogin('jo@gmail.com', 'password')

# Query the server for an Atom feed containing a list of your documents.
#documents_feed = client.GetDocumentListFeed()
# Loop through the feed and extract each document entry.
#for document_entry in documents_feed.entry:
  # Display the title of the docume
#    print document_entry.title.text

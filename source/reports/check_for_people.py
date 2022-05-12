#!/usr/bin/env python3

import json
import csv

people_input_file = "../../config/people.csv"
clean_data_input_file = "../../data/arcnwc/summary_added_to"

# Read in the list of people
people_of_interest = set()
with open(people_input_file, encoding='utf8') as csvfile:
    f = csv.reader(csvfile)
    for row in f:
        people_of_interest.add(row[1] + " " + row[0][0])

print(people_of_interest)

# Read in the data in
papers = []
with open(clean_data_input_file) as jsonfile:
    papers = json.load(jsonfile)
    # print(papers)

# Cycle through each paper checking each one as we go and building some
# html if needed.
html = '<table>'
html += "<tr><th>People of interest</th><th>Full author list</th><th>Title</th><th>DOI</th></tr>"
for this_paper in papers:
    full_author_list = this_paper['clean']['full_author_list']
    # print(full_author_list)
    cleaned_authors = []
    for this_author in full_author_list:
        # print(this_author['clean'])
        cleaned_authors.append(this_author['clean'])

    temp = (people_of_interest.intersection(set(cleaned_authors)))
    # print(temp)
    if len(temp):
        print("Match found")
        print(cleaned_authors)
        print(temp)

        html += "<tr>"

        html += "<td>"
        html += str(temp)
        html += "</td>"

        html += "<td>"
        html += str(cleaned_authors)
        html += "</td>"

        html += "<td>"
        html += this_paper['clean']['title']
        html += "</td>"

        html += "<td>"
        try:
            html += '<a href="http://doi.org/' + this_paper['IDs']['DOI'] + '">' + this_paper['IDs']['DOI'] + '</a>'
        except:
            html += ""
        html += "</td>"

        html += "</tr>"

html += "</table>"

css = """
        <style>
        table, th, td {
            border: 1px solid;
        }
        </style>
        """


html_file = open("people_of_interest.html", 'w')
html_file.write(css)
html_file.write(html)

#! /usr/bin/env python
import json


class ErrorLog:

    errorArray = []
    warningArray = []
    index = 0

    # Used for logging errors and warnings related to the pipeline
    def logError(self, txt):
        self.errorArray.append("[ERROR] " + txt)

    def logWarning(self, txt):
        self.warningArray.append("[WARNING] " + txt)

    # Used for logging errors and warnings relating to a paper
    def logErrorPaper(self, txt, this_paper):
        # Create the text for logging an error
        string = ""

        string += "[ERROR] " + txt + " <strong><a onclick=' var obj = document.getElementById(\"error_extra_info_" + str(self.index) + "\"); if (obj.style.display == \"none\") obj.style.display = \"block\"; else obj.style.display = \"none\";'>Details</a></strong>"

        string += "<div id='error_extra_info_" + str(self.index) + "' style='display:none;'>"
        string += "<table>"
        string += "<tr>"
        string += "<th>Hash</th><th>" + this_paper['IDs']['hash'] + "</th>"
        string += "</tr>"
        string += "<tr>"
        try:
            string += "<td>Title</td><td>" + str(this_paper['title']) + "</td>"
        except:
            string += "<td>Title</td><td>Problem with title</td>"
        string += "</tr>"
        string += "<tr>"
        if this_paper['IDs']['DOI'] != "":
            string += "<td>DOI</td><td><a href='http://doi.org/" + str(this_paper['IDs']['DOI']) + "'>" + str(this_paper['IDs']['DOI']) + "</a></td>"
        else:
            string += "<td>DOI</td><td>No DOI</td>"
        string += "</tr>"
        string += "<tr>"
        try:
            string += "<td>PUBMED ID</td><td><a href='https://www.ncbi.nlm.nih.gov/pubmed/" + str(this_paper['IDs']['PMID']) + "'>" + str(this_paper['IDs']['PMID']) + "</a></td>"
        except:
            string += "<td>PUBMED ID</td><td>No PMID</td>"
        string += "</tr>"
        string += "<tr>"
        string += "<td>Object</td><td><code class='prettyprint'><textarea class='textoutput' style='width:100%;min-height:600px;font-size:14px;font-family:\"Courier New\", Courier, monospace'>" + str(json.dumps(this_paper)).replace("<", "&lt;").replace(">", "&gt;") + "</textarea></code></td>"
        string += "</tr>"
        string += "</table>"
        string += "</div>"

        self.errorArray.append(string)
        self.index += 1

    def logWarningPaper(self, txt, this_paper):
        # Create the text for logging a warning

        string = ""

        string += "[WARNING] " + txt + " <strong><a onclick=' var obj = document.getElementById(\"error_extra_info_" + str(self.index) + "\"); if (obj.style.display == \"none\") obj.style.display = \"block\"; else obj.style.display = \"none\";'>Details</a></strong>"

        string += "<div id='error_extra_info_" + str(self.index) + "' style='display:none;'>"
        string += "<table>"
        string += "<tr>"
        string += "<th>Hash</th><th>" + this_paper['IDs']['hash'] + "</th>"
        string += "</tr>"
        string += "<tr>"
        try:
            string += "<td>Title</td><td>" + str(this_paper['title']) + "</td>"
        except:
            string += "<td>Title</td><td>Problem with title</td>"
        string += "</tr>"
        string += "<tr>"
        if this_paper['IDs']['DOI'] != "":
            string += "<td>DOI</td><td><a href='http://doi.org/" + str(this_paper['IDs']['DOI']) + "'>" + str(this_paper['IDs']['DOI']) + "</a></td>"
        else:
            string += "<td>DOI</td><td>No DOI</td>"
        string += "</tr>"
        string += "<tr>"
        try:
            string += "<td>PUBMED ID</td><td><a href='https://www.ncbi.nlm.nih.gov/pubmed/" + str(this_paper['IDs']['PMID']) + "'>" + str(this_paper['IDs']['PMID']) + "</a></td>"
        except:
            string += "<td>PUBMED ID</td><td>No PMID</td>"
        string += "</tr>"
        string += "<tr>"
        string += "<td>Object</td><td><code class='prettyprint'><textarea class='textoutput' style='width:100%;min-height:600px;font-size:14px;font-family:\"Courier New\", Courier, monospace'>" + str(json.dumps(this_paper)).replace("<", "&lt;").replace(">", "&gt;") + "</textarea></code></td>"
        string += "</tr>"
        string += "</table>"
        string += "</div>"

        self.warningArray.append(string)
        self.index += 1

    # Called at the end of the papers.py script. It returns the html error log that is directly written to the page.
    def printLog(self):
        log = ""

        log += "<h2>" + str(len(self.warningArray)) + " Warnings</h2>"
        log += "<div>"
        for warn in self.warningArray:
            log += warn + "<br/>"
        log += "</div>"

        log += "<h2>" + str(len(self.errorArray)) + " Errors</h2>"
        log += "<div>"
        for error in self.errorArray:
            log += error + "<br/>"
        log += "</div>"

        return log

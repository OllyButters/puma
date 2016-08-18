#! /usr/bin/env python


class ErrorLog:

    errorArray = []
    warningArray = []
    index = 0

    # OLD ERROR LOGGING FUNCTIONS - Shouldn't be used anymore
    # def logError(self, txt):
    #   self.errorArray.append("[ERROR] " + txt)

    # def logWarning(self, txt):
    #   self.warningArray.append("[WARNING] " + txt)

    def logErrorPaper(self, txt, this_paper):
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
            string += "<td>DOI</td><td>" + str(this_paper['IDs']['DOI']) + "</td>"
        else:
            string += "<td>DOI</td><td>No DOI</td>"
        string += "</tr>"
        string += "<tr>"
        try:
            string += "<td>PUBMED ID</td><td>" + str(this_paper['PMID']) + "</td>"
        except:
            string += "<td>PUBMED ID</td><td>No PMID</td>"
        string += "</tr>"
        string += "<tr>"
        string += "<td>Object</td><td style='font-size:14px'>" + str(this_paper).replace("<", "&lt;").replace(">", "&gt;") + "</td>"
        string += "</tr>"
        string += "</table>"
        string += "</div>"

        self.errorArray.append(string)
        self.index += 1

    def logWarningPaper(self, txt, this_paper):
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
            string += "<td>DOI</td><td>" + str(this_paper['IDs']['DOI']) + "</td>"
        else:
            string += "<td>DOI</td><td>No DOI</td>"
        string += "</tr>"
        string += "<tr>"
        try:
            string += "<td>PUBMED ID</td><td>" + str(this_paper['PMID']) + "</td>"
        except:
            string += "<td>PUBMED ID</td><td>No PMID</td>"
        string += "</tr>"
        string += "<tr>"
        string += "<td>Object</td><td style='font-size:14px'>" + str(this_paper).replace("<", "&lt;").replace(">", "&gt;") + "</td>"
        string += "</tr>"
        string += "</table>"
        string += "</div>"

        self.warningArray.append(string)
        self.index += 1

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

#! /usr/bin/env python

class ErrorLog:
    
    errorArray = []
    warningArray = []

    def logError(self, txt ):
        self.errorArray.append( "[ERROR] " + txt )

    def logWarning(self, txt ):
        self.warningArray.append( "[WARNING] " + txt )

    def printLog(self):
        log = ""

        log += "<h2>" + str(len(self.warningArray)) + " Warnings</h2>"
        log += "<p>"
        for warn in self.warningArray:
            log += warn + "<br/>"
        log += "</p>"

        log += "<h2>" + str(len(self.errorArray)) + " Errors</h2>"
        log += "<p>"
        for error in self.errorArray:
            log += error + "<br/>"
        log += "</p>"

        return log



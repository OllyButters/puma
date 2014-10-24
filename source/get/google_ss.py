/**
 * Copies the data from each tab into a separte public file
 * Olly Butters - 20/9/14
 */
function exportPmids() {
  
  //Open main spreadsheet with all the data in it
  //var sheet = SpreadsheetApp.getActiveSheet();
  var sheets = SpreadsheetApp.openById('1Zp4w0ZlxN6eK2n0tNBr-CYZNZmxf-9GXkcYs_Terr4E');
  var start_year = 1999;
  var end_year = 2014;
  
  //Create a new spreadsheet named 'pmid'
  var target = SpreadsheetApp.create('pmids');
  
  Logger.log((end_year-start_year)+' years to process');
  
  for (var tab_no=start_year; tab_no <= end_year; tab_no++)
  {
    //Output what we think we are doing.
    Logger.log('Year='+tab_no);

    //Get this sheet based on its name
    var this_sheet = sheets.getSheetByName(tab_no);
            
    //var this_sheet_name = this_sheet.getSheetName();
    SpreadsheetApp.getActiveSpreadsheet().toast('Doing '+tab_no, 'Status', 100);
    
  
    var rows = this_sheet.getDataRange();
    var numRows = rows.getNumRows();
    //var values = rows.getValues();
    
    //row #, col #, # of rows, # of cols
    var values = this_sheet.getSheetValues(2, 2, numRows, 1);
      
    //Append the data from the original one to the new one
    for (var i = 0; i <= numRows - 1; i++) {
      var row = values[i];
      target.appendRow(row);
    }
  }
  
  
  //Give drive a chance to catch up
  SpreadsheetApp.flush();
  Utilities.sleep(2000);
  
  

  //Move the file to the public dir. This is a bit convoluted.
  //var folder = DriveApp.getFoldersByName('public');
  var folder = DriveApp.getFoldersByName('public_data');
  var files = DriveApp.getFilesByName('pmids');  
  var file = files.next();
  var this_folder = folder.next();
  var date_time = Utilities.formatDate(new Date(), "GMT", "yyyy-MM-dd'T'HH:mm:ss");
  var temp = file.makeCopy(date_time, this_folder);
  DriveApp.removeFile(file);
  
  SpreadsheetApp.getActiveSpreadsheet().toast('Finished :)', 'Status', 100);

  
};

/**
 * Add a menu item to the active spreadsheet. Need to reopen spreadsheet to activate 
 *
*/
function onOpen() {
  var spreadsheet = SpreadsheetApp.getActiveSpreadsheet();
  var entries = [{
    name : "Export PMIDs",
    functionName : "exportPmids"
  }];
  spreadsheet.addMenu("Projects DB menu", entries);
};

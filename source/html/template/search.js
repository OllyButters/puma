String.prototype.contains = function(it) { return this.toLowerCase().indexOf(it.toLowerCase()) != -1; };

function search(){
  // This function searches the papers data object for text matching the query.
  // It is called when the user clicks the search button.

  // If the input is blank then do not search because this will return all papers
  if( document.getElementById("search").value == "" ){
    return;
  }

  // Hide results and display message to show that the search is being carried out.
  document.getElementById("searching").style.display="block";
  document.getElementById("search_results").style.display="none";
  document.getElementById("num_search_results").innerHTML = "";

  // Get the papers and exec data
  //var raw_data = document.getElementById("search_data").innerHTML;
  //var data = JSON.parse(raw_data);
  var data = papers;
//  var raw_exec = document.getElementById("exec_list").innerHTML;
//  var exec_list = JSON.parse(raw_exec);

  // Get the input search query
  var query = document.getElementById("search").value;
  var results = "";
  var num_results = 0;

  // Split the search query up into an array of the words
  var query_components = query.split(" ");

  // Check each paper object
  for( i = 0 ; i < data.length; i++ ){

    // Count the number of matched elements
    // If the number of matched elements equals the number of query_components
    // then the paper is a match
    var components_match = 0

    for( j = 0 ; j < query_components.length; j++ ){
      var match = false;

      // Check title
      if( data[i].merged.title.contains( query_components[j] ) ){
        match = true;
      }

      // Check abstract
      try{
        //if(data[i].merged.MedlineCitation.Article.Abstract.AbstractText[0].contains(query_components[j])){
        if(data[i].clean.abstract.contains(query_components[j])){
          match = true;
        }
      } catch(err){}

      // Check subject text
      //try {
      //  for( n = 0; n < data[i].merged.subject.length ; n++ ){
          //if( data[i].merged.subject[n].contains( query_components[j] ) ){
        //  if( data[i].merged.subject[n].contains( query_components[j] ) ){
          //    match = true;
          //}
        //}
      //} catch (err){}

      // Check keywords
      try {
//        for( n = 0; n < data[i].merged.MedlineCitation.MeshHeadingList.length ; n++ ){
  //        if( data[i].merged.MedlineCitation.MeshHeadingList[n].DescriptorName.contains( query_components[j] ) ){
          for( n = 0; n < data[i].clean.mesh.length ; n++ ){
            if( data[i].clean.mesh[n].term.contains( query_components[j] ) ){

            match = true;
          }
        }
      } catch (err){}

      // Check author names
      try {
        for( n = 0; n < data[i].clean.full_author_list.length ; n++ ){
          if( data[i].clean.full_author_list[n].family.contains( query_components[j] ) || data[i].clean.full_author_list[n].given.contains( query_components[j] ) ){
            match = true;
          }
        }
      } catch (err){}

      if( match ){
        components_match += 1;
      }
    }

    if( components_match == query_components.length ){
      // This is a match write the paper to page

      results += "<div class='paper'>";

        // altmetric data
        try{
            if ( data[i].IDs.DOI ) {
              results += '<div style="float:right;" data-badge-popover="right" data-badge-type="donut" data-doi="' + this_paper['IDs']['DOI'] + '" data-hide-no-mentions="true" class="altmetric-embed"></div>';
            }
        } catch (err){}

        // Paper title
        results += '<span style="text-decoration: underline; font-weight:bold;">' + data[i].clean.title + '</span><br/>';

        // Authors
          authors = [];
          for( a = 0 ; a < data[i].clean.full_author_list.length; a++ ){

              if( a > 0 ){
                results += "; ";
              }
              results += data[i].clean.full_author_list[a].family + ', ' + data[i].clean.full_author_list[a].given;
          }
          results += '<br/>'


        // Journal, volume and issue
        try{
            //if( "ISOAbbreviation" in data[i].merged.MedlineCitation.Article.Journal ){
              results += data[i].clean.journal.journal_name;
            //}
        } catch (err){}

        try{
            if( "volume" in data[i].clean.journal ){
              results += ', Volume ' + data[i].clean.journal.volume;
            }
        } catch (err){}

        try{
            if( "issue" in data[i].clean.journal){
              results += ', Issue ' + data[i].clean.journal.issue;
            }
        } catch (err){}

        results += '<br/>';

        // PMID
        try{
            if (data[i].IDs.PMID){
              results += 'PMID: <a href="http://www.ncbi.nlm.nih.gov/pubmed/' + data[i].IDs.PMID + '">' + data[i].IDs.PMID + '</a>&nbsp;';
            }
        } catch (err){}

        // DOI
        try{
            if( data[i].IDs.DOI ){
                results += 'DOI: <a href="https://doi.org/' + data[i].IDs.DOI + '">' + data[i].IDs.DOI + '</a>&nbsp;';
            }
        } catch (err){}

        // Citation Counts and Sources
        number_citations_counts = 2; // The number of different citation count sources
        citations_counts_width = 100 / number_citations_counts;
        results += "<table class='citation_table'>";
        results += '<tr><th colspan="' + number_citations_counts + '">Citation Counts</th></tr>';
        results += '<tr>';
        try{
            // Try to display citation count with link to scopus page
            results += '<td style="width:' + citations_counts_width + '%;">Scopus: <a href="https://www.scopus.com/record/display.uri?eid=' + data[i].IDs.scopus + '&origin=inward&txGid=0">' + data[i].clean.citations.scopus.count + '</a></td>';
        } catch (err){
            try {
                results += '<td style="width:' + citations_counts_width + '%;">Scopus: ' + data[i].clean.citations.scopus.count + '</td>';
            } catch (err){
                results += '<td style="width:' + citations_counts_width + '%;">Scopus: -</td>';
            }
        }

        try{
            results += '<td style="width:' + citations_counts_width + '%;">Europe PMC: ' + data[i].clean.citation.PMC.count + '</td>';
        } catch (err){
            results += '<td style="width:' + citations_counts_width + '%;">Europe PMC: -</td>';
        }

        results += '</tr>';
        results += "</table>";
        results += '</div>';


      results += "</div>";
      num_results++;
    }

  }

  // Write results to screen and display them and the papers matching count. Also hide the searching message.
  document.getElementById("search_results").innerHTML = results;
  document.getElementById("searching").style.display="none";
  document.getElementById("search_results").style.display="block";
  document.getElementById("num_search_results").innerHTML = "Matches " + num_results + " of " + data.length + " Papers";

}

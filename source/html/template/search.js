String.prototype.contains = function(it) { return this.toLowerCase().indexOf(it.toLowerCase()) != -1; };

function search(){

  document.getElementById("searching").style.display="block";
  document.getElementById("search_results").style.display="none";
  document.getElementById("num_search_results").innerHTML = "";

  var raw_data = document.getElementById("search_data").innerHTML;
  var data = JSON.parse(raw_data);
  var raw_exec = document.getElementById("exec_list").innerHTML;
  var exec_list = JSON.parse(raw_exec);

  var query = document.getElementById("search").value;
  var results = "";
  var num_results = 0;

  var query_components = query.split(" ");

  for( i = 0 ; i < data.length; i++ ){

    var components_match = 0
    for( j = 0 ; j < query_components.length; j++ ){
      var match = false;

      if( data[i].title.contains( query_components[j] ) ){
        match = true;
      }

      try{
        if(data[i].MedlineCitation.Article.Abstract.AbstractText[0].contains(query_components[j])){
          match = true;
        }
      } catch(err){}

      try {
        for( n = 0; n < data[i].subject.length ; n++ ){
          if( data[i].subject[n].contains( query_components[j] ) ){
            match = true;
          }
        }
      } catch (err){}

      try {
        for( n = 0; n < data[i].MedlineCitation.MeshHeadingList.length ; n++ ){
          if( data[i].MedlineCitation.MeshHeadingList[n].DescriptorName.contains( query_components[j] ) ){
            match = true;
          }
        }
      } catch (err){}

      try {
        for( n = 0; n < data[i].author.length ; n++ ){
          if( data[i].author[n].family.contains( query_components[j] ) || data[i].author[n].given.contains( query_components[j] ) ){
            match = true;
          }
        }
      } catch (err){}

      if( match ){
        components_match += 1;
      }
    }

    if( components_match == query_components.length ){
      results += "<div class='paper'>";

      // altmetric data
        try{
            if ( data[i].IDs.DOI ) {
              results += '<div style="float:right;" data-badge-popover="right" data-badge-type="donut" data-doi="' + this_paper['IDs']['DOI'] + '" data-hide-no-mentions="true" class="altmetric-embed"></div>';
            }
        } catch (err){}

        // Paper title
        results += '<span style="text-decoration: underline; font-weight:bold;">' + data[i].title + '</span><br/>';

        // Authors
          authors = [];
          author_on_exec = false;
          for( a = 0 ; a < data[i].author.length; a++ ){
            //Check if an author was on the exec Comittee
              for ( x = 0; x < exec_list.length; x++ ){
                //Check if authors name matches
                if (exec_list[x][2] == data[i].author[a].clean){
                  // Get start and end date of exec membership
                  var start_split = exec_list[x][0].split("/");
                  var end_split = exec_list[x][1].split("/");

                  var exec_start = new Date( start_split[2], start_split[1], start_split[0] ).getTime()/1000;
                  var exec_end = Date.now()/1000;
                  if( exec_list[x][1] != ""){
                    exec_end = new Date(end_split[2], end_split[1], end_split[0]).getTime()/1000;
                  }

                  // Convert issued date into a timestamp
                  var clean_date = data[i].Extras.CleanDate;
                  var issued_timestamp = new Date(clean_date.year,clean_date.month,clean_date.day).getTime()/1000;

                  // If publication is issued between exec_start and exec_end then flag
                  if (issued_timestamp > exec_start && issued_timestamp < exec_end){
                    author_on_exec = true;
                    break;
                  }
                }
              }

              if( a > 0 ){
                results += "; ";
              }
              results += data[i].author[a].family + ', ' + data[i].author[a].given;
          }
          results += '<br/>'

          if( author_on_exec ){
            results += '<div style="text-align:center;font-size:14px;background:#' + secondary_colour + ';color:#' + primary_colour + ';padding:2px 4px;box-shadow: 0px 0px 1px #4e4e4e inset;">At least one author was a member of the ' + name + ' Executive Committee.</div>';
          }


        // Journal, volume and issue
        try{
            if( "ISOAbbreviation" in data[i].MedlineCitation.Article.Journal ){
              results += data[i].MedlineCitation.Article.Journal.ISOAbbreviation;
            }
        } catch (err){}

        try{
            if( "volume" in data[i] ){
              results += ', Volume ' + data[i].volume;
            }
        } catch (err){}

        try{
            if( "Issue" in data[i].MedlineCitation.Article.Journal.JournalIssue ){
              results += ', Issue ' + data[i].MedlineCitation.Article.Journal.JournalIssue.Issue;
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
                results += 'DOI: <a href="http://doi.org/' + data[i].IDs.DOI + '">' + data[i].IDs.DOI + '</a>&nbsp;';
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
            results += '<td style="width:' + citations_counts_width + '%;">Scopus: <a href="https://www.scopus.com/record/display.uri?eid=' + data[i].Extras.eid + '&origin=inward&txGid=0">' + data[i].Extras.Citations + '</a></td>';
        } catch (err){
            try {
                results += '<td style="width:' + citations_counts_width + '%;">Scopus: ' + data[i].Extras.Citations + '</td>';
            } catch (err){
                results += '<td style="width:' + citations_counts_width + '%;">Scopus: -</td>';
            }
        }

        try{
            results += '<td style="width:' + citations_counts_width + '%;">Europe PMC: ' + data[i]["Extras"]["Citations-EuropePMC"] + '</td>';
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

  document.getElementById("search_results").innerHTML = results;
  document.getElementById("searching").style.display="none";
  document.getElementById("search_results").style.display="block";
  document.getElementById("num_search_results").innerHTML = "Matches " + num_results + " of " + data.length + " Papers";

}

google.load("visualization", "1", {packages:["corechart"]});
      google.setOnLoadCallback(drawChart);
      function drawChart() {

        //number per year
        var data = google.visualization.arrayToDataTable(papers);
        var options = {
          title: 'Number of papers published each year for this ' + title,
          colors: [primary_colour],
          vAxis: {minValue: 0,
                 format:"#"},
          legend: {position: 'none'}
        };
        var chart = new google.visualization.ColumnChart(document.getElementById('papers_chart_div'));
        chart.draw(data, options);


 	//number per year
        var data = google.visualization.arrayToDataTable(citations);
        var options = {
          title: 'Number of citations of publications in year for this ' + title,
          colors: [secondary_colour],
          vAxis: {minValue: 0,
                  format:"#"},
          legend: {position: 'none'}
        };
        var chart = new google.visualization.ColumnChart(document.getElementById('citations_chart_div'));
        chart.draw(data, options);


      }


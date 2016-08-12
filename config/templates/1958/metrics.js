google.load("visualization", "1", {packages:["corechart"]});
      google.setOnLoadCallback(drawChart);
      function drawChart() {

        //cumulative
        var data = google.visualization.arrayToDataTable(cumulative);
        var options = {
          title: 'Cumulative number of papers published',
          colors: ["#1d294d"]
        };
        var chart = new google.visualization.LineChart(document.getElementById('cumulative_div'));
        chart.draw(data, options);

        //number per year
        var data = google.visualization.arrayToDataTable(papers_per_year);
        var options = {
          title: 'Number of papers published',
          colors: ["#1d294d"]

        };
        var chart = new google.visualization.ColumnChart(document.getElementById('papers_per_year_div'));
        chart.draw(data, options);

	//Number of Papers for Citation Counts
        var data = google.visualization.arrayToDataTable(papers_per_citation_count);
        var options = {
          title: 'Number of Papers for Citation Counts',
          colors: ["#1d294d"],
	  hAxis: {title:"Number of Citations"}

        };
        var chart = new google.visualization.ColumnChart(document.getElementById('papers_per_citation_count_div'));
        chart.draw(data, options);


      }

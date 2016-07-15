google.load("visualization", "1", {packages:["corechart"]});
      google.setOnLoadCallback(drawChart);
      function drawChart() {

        //cumulative
        var data = google.visualization.arrayToDataTable(cumulative);
        var options = {
          title: 'Cumulative number of papers published',
          colors: ["#c9002f"]
        };
        var chart = new google.visualization.LineChart(document.getElementById('cumulative_div'));
        chart.draw(data, options);

        //number per year
        var data = google.visualization.arrayToDataTable(papers_per_year);
        var options = {
          title: 'Number of papers published',
          colors: ["#c9002f"]

        };
        var chart = new google.visualization.ColumnChart(document.getElementById('papers_per_year_div'));
        chart.draw(data, options);


      }

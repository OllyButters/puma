google.load("visualization", "1", {packages:["corechart"]});
      google.setOnLoadCallback(drawChart);
      function drawChart() {

        //number per year
        var data = google.visualization.arrayToDataTable(papers);
        var options = {
          title: 'Number of papers published each year for keyword',
          colors: ["#c9002f"]

        };
        var chart = new google.visualization.ColumnChart(document.getElementById('papers_chart_div'));
        chart.draw(data, options);


 	//number per year
        var data = google.visualization.arrayToDataTable(citations);
        var options = {
          title: 'Number of citations of publications in year for keyword',
          colors: ["#FFB612"]

        };
        var chart = new google.visualization.ColumnChart(document.getElementById('citations_chart_div'));
        chart.draw(data, options);


      }


    function initialize() {
      var latlng = new google.maps.LatLng(26, 24);
      var myOptions = {
        zoom: 3,
        center: latlng,
        mapTypeId: google.maps.MapTypeId.ROADMAP
      };
      var map = new google.maps.Map(document.getElementById("map_canvas"), myOptions);

      for (var i = 0; i < locations.length; i++) {
  	var site = locations[i];
  	var myLatLng = new google.maps.LatLng(site[1], site[2]);
          var marker = new google.maps.Marker({
              position: myLatLng,
              map: map,
              title: site[0]
          });
      }
  }

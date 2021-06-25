function update_timestamp(timestamp_path) {
    fetch(timestamp_path)
    .then(response => {
        return response.text()
    })
    .then(data => {
        document.getElementById("update_timestamp").innerHTML = "Updated: " + data;
    });
}
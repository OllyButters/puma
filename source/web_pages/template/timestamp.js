function update_timestamp() {
    fetch("https://ollybutters.github.io/puma/timestamp.html")
    .then(response => {
        return response.text()
    })
    .then(data => {
        document.getElementById("update_timestamp").innerHTML = "Updated: " + data;
    });
}
function update_timestamp(timestamp_path) {
    const myInit = {
        mode: 'no-cors'
    };
    fetch(timestamp_path, myInit)
    .then(response => {
        return response.text()
    })
    .then(data => {
        document.getElementById("update_timestamp").innerHTML = "Updated: " + data;
    });
}


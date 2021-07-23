function update_timestamp(timestamp_path) {
    fetch(timestamp_path)
    .then(response => {
        return response.text()
    })
    .then(data => {
        document.getElementById("update_timestamp").innerHTML = "Updated: " + data;
    });
}

// Add a listener to post to parent window the size of the screen. This is used
// when displayed in an iframe and needs the dimensions to resize
window.addEventListener('load', function() {
    //let message = { height: document.body.scrollHeight, width: document.body.scrollWidth };	
    
    content_height = document.getElementById("wrapper").scrollHeight;
    footer_height = document.getElementsByClassName("foot")[0].scrollHeight;

    let message = { height: Number(content_height) + Number(footer_height) + Number(100), width: document.body.scrollWidth };	

	// window.top refers to parent window
	window.top.postMessage(message, "*");
    console.log(message)
});

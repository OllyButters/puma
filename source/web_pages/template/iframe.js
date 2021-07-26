if ( window.location !== window.parent.location )
{
    
    // The page is in an iFrames
    console.log("The page is in an iFrame");
    
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

        // hide the header
        var elem = document.querySelector('#header-container');
        elem.style.display = 'none';

    });

} 
else {
        
    // The page is not in an iFrame
    console.log("The page is not in an iFrame");
}





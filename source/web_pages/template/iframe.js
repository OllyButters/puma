if ( window.location !== window.parent.location )
{
    // The page is in an iFrames
    console.log("The page is in an iFrame");
    
    function resize_parent() {

        // hide the header
        var elem = document.querySelector('#header-container');
        elem.style.display = 'none';

        content_height = document.getElementById("wrapper").scrollHeight;
        footer_height = document.getElementsByClassName("foot")[0].scrollHeight;

        let message = { height: Number(content_height) + Number(footer_height) + Number(100), width: document.body.scrollWidth };	

        // window.top refers to parent window
        window.top.postMessage(message, "*");
        console.log(message);

    }

    // Add a listener to post to parent window the size of the screen. This is used
    // when displayed in an iframe and needs the dimensions to resize
    window.addEventListener('load', function() {

        // Call the resize_parent a few times with a gap between them.
        // It's an ugly hack as theres loads of async calls to external
        // services which I don't have call backs for. Sorry.
        resize_parent();
        setTimeout( function(){ resize_parent(); }, 1000); 
        setTimeout( function(){ resize_parent(); }, 5000); 
    });

} 
else {
    // The page is not in an iFrame
    console.log("The page is not in an iFrame");
}


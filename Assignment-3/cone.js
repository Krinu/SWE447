var gl = null;

function init() {
    var canvas = document.getElementById( "webgl-canvas" );

    gl = WebGLUtils.setupWebGL( canvas );

    if ( !gl ) {
        alert("Unable to setup WebGL");
        return;
    }

    gl.clearColor( 0.0, 0.3, 0.3, 1.0 );

    render();
}

function render() {
    gl.clear( gl.COLOR_BUFFER_BIT );
}

window.onload = init;

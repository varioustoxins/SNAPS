$(function () {
    insertDownloadLinks(files);
});

function insertDownloadLinks(files) {
    if (files) {
        $("#files").append("<p>Latest run:</p>");

        if (files['results'])
            $("#files").append("<button id='downloadResults'>Download results</button>");
        if (files['plot'])
            $("#files").append("<button id='downloadPlot'>Download plot</button>");
    }
    $("#downloadResults").click(function () {
        download(files['results']);
    });
    $("#downloadPlot").click(function () {
        download(files['plot']);
    });
}

function download(filePath) {
    $.ajax({
        type: "POST",
        url: $SCRIPT_ROOT + '/download',
        data: { filePath: filePath },
        success: function (output) {
            var naps_download = new Blob([output], { type: "application / octet - stream" });
            saveAs(naps_download, filePath.split("\\")[1]);
        }
    });
}
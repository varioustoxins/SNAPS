function insertDownloadLinks(files) {
    $("#downloadResults").empty();
    $("#downloadPlot").empty();
//    $("#email").empty();

    if (!jQuery.isEmptyObject(files)) {
        if (files['results'])
//            $("#downloadResults").append("<img src='/static/save.png' id='downloadResultsButton' class='fileDownload' title='Download results'/>");
            $("#downloadResults").append("<button id='downloadResultsButton'>Download results table</button>");
        if (files['plot'])
//            $("#downloadPlot").append("<img src='/static/save.png' id='downloadPlotButton' class='fileDownload' title='Download plot'/>");
            $("#downloadPlot").append("<button id='downloadPlotButton'>Download results strip plot</button>");
        if (files['log'])
            $("#downloadLog").append("<button id='downloadLogButton'>Download logs</button>");

//        if (files['results'] || files['plot'])
//            $("#email").append("<form id='emailForm'><input type='text' name='emailAddress' placeholder='your@email.com' class='fileSender'><button id='emailSubmit'>Email all results</button></form>");
    }
    $("#downloadResultsButton").click(function () {
        download('results', files['results']);
    });
    $("#downloadPlotButton").click(function () {
        download('plot', files['plot']);
    });
    $("#downloadLogButton").click(function () {
        download('log', files['log']);
    });
/*    $(document).on("submit", "#emailForm", function (event) {
        event.preventDefault();
        var data = new FormData(this);
        $.each(files, function (key, value) {
            if (value)
                data.append(key, value);
        });
        $.ajax({
            url: $SCRIPT_ROOT + '/email',
            type: "POST",
            dataType: "JSON",
            data: data,
            processData: false,
            contentType: false,
            success: function (data) {
                if (data.status !== 'ok') {
                    alert(data.message);
                }
                else {
                    alert("Results emailed successfully");
                }
            },
            error: function (err) {
                alert("Email failed");
                console.log(err);
            }
        });
    });*/
}

function download(fileName, file) {
    var naps_download = new Blob([file], { type: "application / octet - stream" });
    saveAs(naps_download, fileName);
}
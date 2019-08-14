function insertDownloadLinks(files) {
    $("#downloadResults").empty();
//    $("#downloadStripPlot").empty();

    if (!jQuery.isEmptyObject(files)) {
        if (files['results'])
            $("#downloadResults").append("<button id='downloadResultsButton'>Download results table</button><br>");
        if (files['shiftlist'])
            $("#downloadResults").append("<button id='downloadShiftlistButton'>Download assigned chemical shifts</button><br>");
        if (files['hsqc_plot'])
            $("#downloadResults").append("<button id='downloadHsqcPlotButton'>Download HSQC plot</button><br>");
        if (files['strip_plot'])
            $("#downloadResults").append("<button id='downloadStripPlotButton'>Download strip plot</button><br>");
        if (files['log_file'])
            $("#downloadResults").append("<button id='downloadLogButton'>Download log file</button><br>");

    }
    $("#downloadResultsButton").click(function () {
        download('results', files['results']);
    });
    $("#downloadShiftlistButton").click(function () {
        download('shiftlist', files['shiftlist']);
    });
    $("#downloadHsqcPlotButton").click(function () {
        download('hsqcPlot.htm', files['hsqc_plot_file']);
    });
    $("#downloadStripPlotButton").click(function () {
        download('stripPlot.htm', files['strip_plot_file']);
    });
    $("#downloadLogButton").click(function () {
        download('log', files['log_file']);
    });
}

function download(fileName, file) {
    var snaps_download = new Blob([file], { type: "application / octet - stream" });
    saveAs(snaps_download, fileName);
}
$(function () {
    $("#run").click(function () {
        $("#form").append("<div id=runLoading class=\"spinner-border\" role=\"status\"><span class=\"sr-only\">Loading...</span></div>");
        $('#resultsSection').collapse('hide');
        $("#tableData").empty();
        $("#errors").empty();
        $("#stripPlot").empty();
        $("#files").empty();
        $("#log").empty();
        $("#form").submit();
        return false;
    });
    $(document).on("submit", "#form", function (event) {
        event.preventDefault();
        $.ajax({
            url: $SCRIPT_ROOT + '/run',
            type: $(this).attr("method"),
            dataType: "JSON",
            data: new FormData(this),
            processData: false,
            contentType: false,
            success: function (data) {
                success(data);
            },
            error: function (err) {
                console.log(err);
            }
        });
    });
});

function success(data) {
    $("#runLoading").remove();
    if (data.status === 'ok') {
        $("#tableData").empty();
        $("#errors").empty();
        $("#stripPlotTopLevel").replaceWith('<div id="stripPlotTopLevel"><div id="stripPlot"></div></div>');
        $("#files").empty();
        $("#log").empty();
        $("#downloadResults").empty();
        $("#tableTopLevel").replaceWith("<div id='tableTopLevel'><table id='table' data-height='500'><thead><tr id='tableData'></tr></thead></table></div>");
        $('#resultsSection').collapse('show');
        insertDownloadLinks(data.files);
        
        $("#downloadResults").prepend("<a href='#tableTopLevel'>Jump to assignment results table</a><br>")
        $("#tableTopLevel").prepend("<h4>Assignment results table</h4>")
        $("#tableTopLevel").append("<a href='#resultsSection'>Return to top of results</a>")
        $.each(data.headers, function (index, header) {
            $("#tableData").append("<th data-field=" + header + ">" + header + "</th>");
        });
        $('#table').bootstrapTable({
            data: data.result
        });
        if (data.files.stripPlot) {
            Bokeh.embed.embed_item(data.files.stripPlot, "stripPlot");
            $("#downloadResults").prepend("<a href='#stripPlotTopLevel'>Jump to assignment strip plot</a> | ")
            $("#stripPlotTopLevel").prepend("<h4>Strip plot</h4>")
            $("#stripPlotTopLevel").append("<a href='#resultsSection'>Return to top of results</a>")
        }
    }
    else if (data.status === 'validation_failed') {
        $.each(data.errors, function (index, error) {
            $("#errors").append("<p>" + error + "</p>");
        });
    }
}
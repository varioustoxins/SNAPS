$(function () {
    $("#run").click(function () {
        $("#form").append("<div id=runLoading class=\"spinner-border\" role=\"status\"><span class=\"sr-only\">Loading...</span></div>");
        $('#resultsSection').collapse('hide');
        $("#tableData").empty();
        $("#errors").empty();
        $("#plot").empty();
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
        $('#resultsSection').collapse('show');
        insertDownloadLinks(data.files);
        if (data.files.plot) {
            Bokeh.embed.embed_item(data.files.plot, "plot");
            $("#downloadResults").before("<a href='#plot'>Jump to assignment strip plot</a> | ")
            $("#plot").prepend("<h4>Strip plot</h4>")
            $("#plot").after("<a href='#resultsSection'>Return to top of results</a>")
        }
        $("#downloadResults").before("<a href='#tableTopLevel'>Jump to assignment results table</a>")
        $("#tableTopLevel").prepend("<h4>Assignment results table</h4>")
        $("#tableTopLevel").append("<a href='#resultsSection'>Return to top of results</a>")
        $.each(data.headers, function (index, header) {
            $("#tableData").append("<th data-field=" + header + ">" + header + "</th>");
        });
        $('#table').bootstrapTable({
            data: data.result
        });
    }
    else if (data.status === 'validation_failed') {
        $.each(data.errors, function (index, error) {
            $("#errors").append("<p>" + error + "</p>");
        });
    }
}
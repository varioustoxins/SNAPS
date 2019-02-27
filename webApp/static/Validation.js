$(function () {
    $("#observedShiftsFile").change(function () {
        //validateInput(this.files[0]);
    });
});

function validateInput(file) {
    var reader = new FileReader();
    reader.onload = function (ev) {
        var str = ev.target.result;
    };
    reader.readAsText(file);
}
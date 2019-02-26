$(function () {
    $(".dropdown-menu span").click(function () {
        $(this).parents(".dropdown").find('.btn').html($(this).text());
        $(this).parents(".dropdown").find('input').val($(this).text());
    });
    $("#advanced").click(function () {
        $('#advancedToggle').collapse('toggle');
        return false;
    });
});
var scrollFunc = function (id) {
    var y = $(id).offset().top - 10;
    window.scroll(0, y);
}
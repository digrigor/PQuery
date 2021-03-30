function alert_message() {
  alert("hey alert");
}


$(document).ready(function() {
                    $('#samples-drop-sa').multiselect({
                        includeSelectAllOption: true,
                        includeResetOption: true,
                        enableFiltering: true,
                    });
                });

$(document).ready(function() {
                    $('#samples-drop-co').multiselect({
                        enableFiltering: true
                    });
                });

$(document).ready(function() {
                    $('#genes-drop-gr').multiselect({
                        enableFiltering: true
                    });
                });

$(document).ready(function() {
                    $('#genes-drop-gr2').multiselect({
                        enableFiltering: true
                    });
                });

$(document).ready(function() {
                    $("#shipping_selector").change(function() {
                      var id = $(this).val();
                      $("#data_" + id).show();
                    }).change();
                });


$(".custom-file-input").on("change", function() {
  var fileName = $(this).val().split("\\").pop();
  $(this).siblings(".custom-file-label").addClass("selected").html(fileName);
});

$(function() {

  // We can attach the `fileselect` event to all file inputs on the page
  $(document).on('change', ':file', function() {
    var input = $(this),
        numFiles = input.get(0).files ? input.get(0).files.length : 1,
        label = input.val().replace(/\\/g, '/').replace(/.*\//, '');
    input.trigger('fileselect', [numFiles, label]);
  });

  // We can watch for our custom `fileselect` event like this
  $(document).ready( function() {
      $(':file').on('fileselect', function(event, numFiles, label) {

          var input = $(this).parents('.input-group').find(':text'),
              log = numFiles > 1 ? numFiles + ' files selected' : label;

          if( input.length ) {
              input.val(log);
          } else {
              if( log ) alert(log);
          }

      });
  });

});



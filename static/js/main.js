
$(document).ready(function (e) {   
    
    // Hide the Spinner when the page is loaded
    $('#spinner').hide();
    // Automatically disable the submit button to prevent any accidental submit without a file selected
    $('#upload-button').attr("disabled", true);

 

    // -----------------------------------------------------------
    // File selection and validation based on the extension
    // -----------------------------------------------------------
    $('#input_file').change(function (e) {        
        // Validation of the type of input: making sure that the extension of the file is either csv or tsv
        const allowed_extensions = ["csv", 'tsv'];
        const errror_message_file_type = "Please provide a file with the following extensions: " + allowed_extensions.join(", ");
        const allowed_size = 2000000 // bites
        const errror_message_file_size = "The size of the file is too large for this Demo website, please provide a file with maximum size of 2MB";
        // Get the name of the file and its extension
        var file = e.target.files[0];
        ext = file.name.split(".");
        ext = ext[ext.length - 1];

        // define empty errors array
        var errors = [];

        $("#messages").empty(); 

        // Check if the filetype is allowed based on the file's extension
        if (!allowed_extensions.includes(ext.toLowerCase())) {
            errors.push(errror_message_file_type);            
        }
        // Check if the size of the file
        if (file.size > allowed_size){
            errors.push(errror_message_file_size);       
        } 
        // Check errors array if empty
        if(errors.length > 0){
            $(this).css('border', 'solid 1px red')
            errors.forEach(function(item){
                $("#messages").append('<li class="text-danger">' + item + '</li>');
                
                                  
            });
            $('#upload-button').attr('disabled', true);
        }else {
            $(this).css('border', 'solid 1px green')
            $('#upload-button').attr('disabled', false);
        }
    });

    // -----------------------------------------------------------
    // File upload using Ajax
    // -----------------------------------------------------------
    $("#form").submit(function (e) {
        // Empty all the results in case of submission of a new file
        empty_result_divs();

        e.preventDefault();       
        $.ajax({
            url: $(this).attr('action'), // Get the url of the post request
            // url: '/file_upload/',
            type: "POST",
            data: new FormData(this), // Construct a form object
            contentType: false,
            cache: false,
            processData: false,
            headers: { "X-CSRFToken": $('meta[name=csrf-token]').attr('content') },
            beforeSend: function () {
                $("#spinner").show();
                $("#messages").empty();
            },
            success: function (data) {                

                if (data['status'] == 'success') {
                    // If the submission status is success (see returned object from back-end)                    
                    $.each(data['messages'], function () {
                        $("#messages").append('<li class="text-success">' + data['messages'] + '</li>');
                    });  
                    
                    // Populate the samples drop-down
                    $('#samples-drop-down').empty();
                    $('#samples-drop-down').append(data_to_dropdown(data['samples'], 'Samples', 12, ["form-control", 'form-control-sm', 'select'], 'select-sample'));                    

                    // Clear the input file field and disable the submit button
                    $('#input_file').val('');
                    $('#upload-button').attr("disabled", true);

                } else {
                    
                    // If the submission status is not success (see returned object from back-end)  
                    data['messages'].forEach(function(item){
                        $("#messages").append('<li class="text-danger">' + item + '</li>');
                    })      
                    
                    // Clear the input file field and disable the submit button
                    $('#samples-drop-down').empty();
                    $('#input_file').val('');
                    $('#upload-button').attr("disabled", true);
                }
            },
            error: function (e) {
                // If any other non managed error
                $("#messages").html(e).fadeIn();
            },
            complete: function () {
                $("#spinner").hide();
            }
        }); //ajax call
    });

    
    // ------------------------------------------------------------------------
    // Change the values of taxonomy dropdown depending on the selected taxon
    // ------------------------------------------------------------------------
    $(document).on('change', "#taxonomy-drop-down select.taxonomy", function () {       

        var tax = $(this).val();

        if (tax !== '') {
            change_taxonomy_dropdown_depending_on_taxon(tax);
        };

        // Change the color of selected taxa in the dropdowns
        color_taxons(tax, $(this).attr('id'));

        // Fetch the taxonomy
        var el = $("#samples-drop-down  #select-sample")
        fetch_taxon_by_sample(el.val());  
        taxonomy_breadcrumbs(tax);

        $('#taxonomy-drop-down select.taxonomy').val("");
    });    
    
    // -----------------------------------------------------------
    // Get taxonomy based on the selected sample
    // -----------------------------------------------------------    
    $(document).on('change', "#samples-drop-down", "#select-sample", function () {
       
        // Clear Chart and breadcrumbs and details of taxonomy if any        
        clear_chart();
        $('#taxonomy-breadcrumbs').empty();
        $('#additional_information').empty();

        var el = $("#samples-drop-down  #select-sample")
        
        url = 'taxonomy_by_sample_name/' + el.val();        
        $.ajax({
            url: url,            
            type: "GET",
            contentType: false,
            cache: false,
            processData: false,
            beforeSend: function () {
                $("#spinner").show();
            },
            success: function (data) {
                if (data['status'] == 'success') { 

                    // Populate the taxonomy dropdowns
                    var levels = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Specie"];
                    $('#taxonomy-drop-down').empty();
                    levels.forEach(function (tax) {                        
                        $('#taxonomy-drop-down').append(data_to_dropdown(
                            data['taxonomy'][tax],
                            tax, 2,
                            ["form-control", 'form-control-sm', "taxonomy", tax.toLowerCase()],
                            tax.toLowerCase()));
                    });

                } else {
                    console.log("An unidentified error occured, please refresh your browser.");
                }
            },
            error: function (e) {
                // If any other non managed error
                $("#messages").html(e).fadeIn();
            },
            complete: function () {
                $("#spinner").hide();
            }
        }); //ajax call

    });

    
    // -----------------------------------------------------------
    // Hide the chart legend on small screens
    // -----------------------------------------------------------
    $(window).resize(function () {
        var width = $('#chart-container').width();
        if (width <= 300) {
            if (myNewChart) {
                myNewChart.options.legend.display = false;
                // myNewChart.options.legend.position = 'top';
                myNewChart.update();
            }
        } else {
            if (myNewChart) {
                myNewChart.options.legend.display = true;
                // myNewChart.options.legend.position = 'right';
                myNewChart.update();
            }
        }
    });


    // -------------------------------------------------
    // click on beradcrumbs updates the graph
    // -------------------------------------------------
    $(document).on('click', "#taxonomy-breadcrumbs .breadcrumb .breadcrumb-item a", function (e) {
        e.preventDefault();
        var tax = $(this).attr('tax');       
        
        if (tax !== '') {
            change_taxonomy_dropdown_depending_on_taxon(tax);
        };

        // color_taxons(tax, $(this).attr('id'));

        // Fetch the taxonomy
        var el = $("#samples-drop-down  #select-sample")
        fetch_taxon_by_sample(el.val());  
        taxonomy_breadcrumbs(tax);
        $('#taxonomy-drop-down select.taxonomy').val("");
    
    });


    // -------------------------------------------
    // Get OTU abundance from the select menu
    // -------------------------------------------
    $(document).on('change', '#otu-abundance', function(){
        
        $('#otu-relative-abundance').empty();
        if($(this).val() != ""){
            var idxOtu = currentTaxon['abundance']['non_zero_otus_numbers'].indexOf(parseInt($(this).val()));        
            var ra = currentTaxon['abundance']['non_zero_otus_abundance_global'][idxOtu];

            $('#otu-relative-abundance').append(
                '<p> Relative abundance: <span class="text-success">' + ra.toExponential(2) +'</span></p>');
        }       
        
        // Clear the histogram
        clear_histogram();
        // Get the index of the otu in the labels
        var idxplot = myHistogram.data.labels.indexOf(parseInt($(this).val()));
        // Create a new colors array and change the color at the index of the otu of interest
        var backgroundColor = [];
        for (i=0;i<myHistogram.data.labels.length;i++){
            if(i == idxplot){
                backgroundColor.push('#000');
            }else{
                backgroundColor.push('rgba(255, 99, 132, 1)');
            }
        }

        // Reset the hitogram
        var dataValues = currentTaxon['abundance']['non_zero_otus_abundance_global'];
        var dataLabels = currentTaxon['abundance']['non_zero_otus_numbers'];
        dataValues = dataValues.map(function(x){
            return -1/Math.log10(x); //Get a better representation of the relative abundance
        });
        draw_histogram(dataLabels, dataValues, backgroundColor);
    });

    

});



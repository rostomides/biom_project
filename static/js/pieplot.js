// -------------------------------------------------------------------------------------
// Initializing global variables that will be updated when charts are created on the fly
var canvas = null;
var ctx = null;
var myNewChart = null; 
var myHistogram = null;
var canvas1 = null;
var ctx1 = null;

// -------------------------------------------------------------------------------------
 


// -------------------------------------
// empty all result dis
// -------------------------------------
function empty_result_divs() {
    var divs = ['samples-drop-down', 'taxonomy-drop-down', 'taxonomy-breadcrumbs', 'additional_information'];
    divs.map(function(id){
        $('#' + id ).empty();
    });

    clear_chart();
}
 
// -------------------------------------------------------------------------------------
// Function setup and draw pie
// -------------------------------------------------------------------------------------
function setup_and_draw_chart(taxonomy, abundances, taxon_level_name, width, height) {
    // Clear the existing chart 
    try {
        clear_chart()
    }
    catch {
        console.log();
    }    

    // Initilize a chart
    canvas = document.getElementById("myChart");
    ctx = canvas.getContext('2d');
    // myNewChart = draw_pie(canvas, ctx, data);
    myNewChart = draw_pie(canvas, ctx, taxonomy, abundances, taxon_level_name)

    // Attribute a click event to the newly created chart
    // This function is applicable if and only if a new canevas is creaed
    canvas.onclick = function (evt) {
        var activePoints = myNewChart.getElementsAtEvent(evt);

        // If the portion defining a taxon is clicked bring the corresponding graph
        if (activePoints[0]) {
            var chartData = activePoints[0]['_chart'].config.data;
            var idx = activePoints[0]['_index'];
            var label = chartData.labels[idx];
            var value = chartData.datasets[0].data[idx];    

            // If a portion of the graph is clicked adjust the grap depensing on the slected taxon
            if(taxon_level_name.toLowerCase() != "specie"){
                // If we get the label from the plot only it will be trucated to the last taxon
                // We will get the last taxon from the breadcrumbs

                var lasteLevel = $("#taxonomy-breadcrumbs .breadcrumb .breadcrumb-item:last-of-type a").attr('tax');
                var taxon = lasteLevel + "; " + label;

                update_chart_after_clicking_portion(taxon);  
                color_taxons(taxon, taxon_level_name);              
            }            
        }
    }; 
};

function update_chart_after_clicking_portion(taxon){   
   
    change_taxonomy_dropdown_depending_on_taxon(taxon);
    // Fetch the taxonomy
    var el = $("#samples-drop-down  #select-sample")
    fetch_taxon_by_sample(el.val());
    taxonomy_breadcrumbs(taxon);
    
    $('#taxonomy-drop-down select.taxonomy').val("");
    
}

// -----------------------------------------------------------
// Populate data into dropdown and append them into a div
// -----------------------------------------------------------
function data_to_dropdown(data, label, numColumnsMd, cssClasses = [], id) {
    var ddown = ""
    ddown += '<div class="form-group col-md-' + numColumnsMd + '">'
    ddown += '<label>' + label + '</label>';
    if (cssClasses.length > 0) {
        ddown += '<select id="' + id + '" class="' + cssClasses.join(" ") + '"><option selected></option>';
    } else {
        ddown += '<select id="' + id + '"><option  selected></option>';
    }

    data.forEach(function (item) {
        ddown += '<option ' + 'value="' + item + '">' + item + '</option>';
    });
    ddown += '</select>';
    ddown += '</div>';
    return ddown;
}

// -----------------------------------------------------------
// Color taxon from lower levels that belong to the selected taxon
// ----------------------------------------------------------- 
function color_taxons(taxon, level){
    var levels = ["kingdom", "phylum", "class", "order", "family", "genus", "specie"];
    var indexOfLevel = levels.indexOf(level.toLowerCase());
    
    $('option').removeClass('tax-selected');

    for(i=0; i<levels.length; i++){

        // For heigher level of taxonomy, the taxons will contan the whole name
        if (i > indexOfLevel){                
            $("#taxonomy-drop-down #" + levels[i] + " option[value^='" + taxon + "']").each(function(){
                $(this).addClass('tax-selected');
            });
        }
        // For lower the will containe just a portion of it
        else {
            var portion = taxon.split("; ").slice(0, i + 1).join("; ");
            $("#taxonomy-drop-down #" + levels[i] + " option[value^='" + portion + "']").each(function(){
                $(this).addClass('tax-selected');
            });
        }
    }
}


// ------------------------------------------------------------------------------------
// Change the content of the dropdowns depending on the last selected taxon
// ------------------------------------------------------------------------------------
function change_taxonomy_dropdown_depending_on_taxon(taxon){
    var levels = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Specie"];
    // the the value selected is not empty
    var taxonomy = taxon.split('; ');
    // Split the selected taxon name
    for (i = 0; i < levels.length; i++) {
        if (i > taxonomy.length - 1) {
            $("#" + levels[i].toLocaleLowerCase()).val("");
        } else {
            // Loop through the levels and each time join the taxonomy array starting from 0 to the index of the taxonomy level
            var tax_pattern = taxonomy.slice(0, i + 1).join("; ");
            // Select the level in the corresponding select drop-down
            // Remove selected from the former selected option                    
            $("#" + levels[i].toLocaleLowerCase() + ' option[value="' + tax_pattern + '"]').prop("selected", true);            
        }
     }
}

// ---------------------------------------------------------------------------------------------------------
// Set the taxonomy in bread crumbs
// ---------------------------------------------------------------------------------------------------------
function taxonomy_breadcrumbs(taxon){
    levels = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Specie"]

    $('#taxonomy-breadcrumbs').empty();
    
    var taxonomy = taxon.split("; ");

    var html = '<nav aria-label="breadcrumb"><ol class="breadcrumb">';

    for(i = 0; i < taxonomy.length -1 ; i++){
        html +=  '<li class="breadcrumb-item"><a href="#" tax="' + taxonomy.slice(0, i+1).join("; ") + '">' + taxonomy[i] + '</a></li>';
    };

    html +=  '<li class="breadcrumb-item"><a class="bc-link" tax="' + taxonomy.slice(0, taxonomy.length).join("; ") + '">' + taxonomy[taxonomy.length -1] + '</a></li></ol></nav>';

    $('#taxonomy-breadcrumbs').append(html);    

}


// ---------------------------------------------------------------------------------------------------------
// Get the data related to the taxon depending on the sample, level of taxonomy and taxon name
// And plot them
// ---------------------------------------------------------------------------------------------------------

// This is a setter for the global variable that will containe the data
function set_global(d){
    currentTaxon = d;
}

function fetch_taxon_by_sample(sample){
    var height = 400;
    var width = 400;

    // Get the lowest taxonomy level selected
    var levels = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Specie"];
    var lastTaxLevel = "";
    levels.forEach(function(level){
        if($('#' + level.toLocaleLowerCase()).val() != ""){
            lastTaxLevel = level.toLocaleLowerCase();
        }
        return;
    });

    // If no taxon selected =>don't do anything
    if (lastTaxLevel ==""){
        return;
    }
    
    // Else proceed with the ajax call
    // Get the name of the taxon and remove spaces from it
    var lastTaxon = $("#" + lastTaxLevel).val().replace(/ /g,'');  
      
    // Perform Ajax call
    var url =  'taxon_abundance_by_sample/' + sample + '/' + lastTaxLevel + "/" + lastTaxon;  
    
    $.ajax({
        url: url,            
        type: "GET",
        contentType: false,
        cache: false,
        processData: false,
        dataType: 'json',            
        beforeSend: function () {
            $("#spinner").show();
        },
        success: function (data) {
            if (data['status'] == 'success') {                
                console.log(data);
                setup_and_draw_chart(
                    data['abundance']['taxonomy_of_next_level'], 
                    data['abundance']['relative_abundances_next_level_global'],                        
                    data['abundance']["main_taxon_level_name"],                         
                    width, height);
                
                // Store the ajax data in the global variables
                set_global(data) ; 
                

                $("#additional_information").empty();
                
                 var otuNumbers = '<select class="form-control form-control-sm" id="otu-abundance" ><option></option>';   
                 data['abundance']["non_zero_otus_numbers"].sort().map(function(item){
                    otuNumbers += '<option>' + item + '</option>';
                 });
                 otuNumbers += '</select>';
                 


                $("#additional_information").append(
                    '<div class="card"> <div class="card-body"> <ul class="list-group list-group-flush"> <li class="list-group-item">  <div> <canvas id="smallChart" width="100px" height="100px" aria-label="blah blah" role="img"></canvas> </div> </li> <li class="list-group-item">  <p>  Number of OTUs : <span class="text-success">' + 
                                                data['abundance']['non_zero_otus_abundance_global'].length + 
                                            '</span> <p>' +  otuNumbers + '<div id="otu-relative-abundance"></div> </li> </ul> </div> </div> ');


                // Draw the histogram
                canvas1 = document.getElementById("smallChart");
                ctx1 = canvas1.getContext('2d');
                var dataValues = data['abundance']['non_zero_otus_abundance_global'];
                var dataLabels = data['abundance']['non_zero_otus_numbers'];

                dataValues = dataValues.map(function(x){
                    return -1/Math.log10(x);
                })

                draw_histogram(dataLabels, dataValues);
                

                                   
            } else {
                console.log("An unidentified error occured, please refresh your browser.");
            }
        },
        error: function (e) {
            // If any other non managed error
            $("#messages").html(e).fadeIn();
        },
        complete: function (data) {
            $("#spinner").hide();
            
        }
    });//ajax call     
}

function draw_histogram(dataLabels, dataValues, backgroundColor='rgba(255, 99, 132, 1)'){
    myHistogram = new Chart(ctx1, {
        type: 'bar',
        data: {
            labels: dataLabels,
            datasets: [{
            // label: 'Group A',
            data: dataValues,
            backgroundColor: backgroundColor ,
            }]
        },
        options: {
            title: {
                display: true,
                text: 'OTU Distribution'
            },
            legend: {
                display: false                            
            },
            scales: {
                xAxes: [{
                    display: false,
                    barPercentage: 1,
                    ticks: {
                        max: 3,
                    },
                    scaleLabel: {
                        display: true,
                        labelString: 'Samples'
                    }
                }, {
                    display: false,                            
                }],
                yAxes: [{
                    ticks: {
                        stepSize: 1,
                        beginAtZero:true
                    },
                    scaleLabel: {
                        display: true,
                        labelString: '-1/log10(Relative Abundance)'
                    }
                }]
            }
        }
    });
}

function clear_histogram() {    
    if(myNewChart){
        ctx1.clearRect(0, 0, canvas1.width, canvas1.height);        
        myHistogram.destroy();
    }
}

// -------------------------------------------------------------------------------------
// Clear the chart 
// -------------------------------------------------------------------------------------   
function clear_chart() {
    // $("#clear-canevas").click(function () {
    if(myNewChart){
        ctx.clearRect(0, 0, canvas.width, canvas.height);
        // myNewChart = null;
        // canvas = null;
        // ctx = null;
        // $("#chart-container").empty();
        myNewChart.destroy();
    }
}

// -------------------------------------------------------------------------------------
// Generate colors for plot
// -------------------------------------------------------------------------------------    
function generate_colors(data) {
    var colorPalette = [];
    for (var i = 0; i < data.length; i++) {
        colorPalette.push("#" + Math.floor(Math.random() * 16777215).toString(16));
    }
    return colorPalette;
};


// -------------------------------------------------------------------------------------
// Draw chart function
// -------------------------------------------------------------------------------------
function draw_pie(canvas, ctx, taxonomy, abundances, taxon_level_name) {

    // Generate background colors
    bgColors = generate_colors(taxonomy);
    borderColors = generate_colors(taxonomy);

    var labels = taxonomy.map(function(item){
        var tax = item.split('; ');
        return tax[ tax.length -1 ];
    });
    var abundances = abundances;

    // Use circles in legend instead of rectangles
    Chart.defaults.global.legend.labels.usePointStyle = true;

    // Seeting the chart by changing the default behaviour of an onClick event in pie chart specifically
    var original = Chart.defaults.pie.legend.onClick;
    Chart.defaults.pie.legend.onClick = function(e, legendItem) {    
    // The label in the chart correspond only to the last porttion of the taxonomy
    // To get the rest target the last part of the bread crumbs
    var lasteLevel = $("#taxonomy-breadcrumbs .breadcrumb .breadcrumb-item:last-of-type a").attr('tax');
    var taxon = lasteLevel + "; " + legendItem.text;

        // Check if OTU is in the label
        if(legendItem.text.indexOf('OTU_') < 0 ){ 
            //  If so update the chart
            update_chart_after_clicking_portion(taxon);  
            color_taxons(taxon, taxon_level_name);
        }else {
            // Prevent click events altogether
            e.stopPropagation();
        }
        // original.call(this, e, legendItem); //Useful in other cases keep it only for documentation
    };

    var myNewChart = new Chart(ctx, {

        plugins: [{
            beforeInit: function(chart, options) {
              chart.legend.afterFit = function() {
                this.height = this.height + 1000;
              };
            }
          }],
        type: 'pie'
        , data: {
            labels: labels,
            datasets: [{
                label: "Abundance of " +  taxon_level_name,
                data: abundances,
                backgroundColor: bgColors,
                borderColor: borderColors,
                borderWidth: 1
            }]
        },
        options: {
            responsive: true, //needed if the chart has to scale to the canvas element
            maintainAspectRatio: false,
            legend: {
                display: true,
                position: 'right',                
            }            
        }
    });

    

    return (myNewChart)
} //function draw_pie(canvas, ctx)








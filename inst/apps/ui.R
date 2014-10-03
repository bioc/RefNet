library(shiny)

# load("tbl.alk.RData")
# providers <- c("ANY", sort(unique(tbl$provider)))

printf <- function(...)print(noquote(sprintf(...)))

providers <-  c("ALL providers",
                "gerstein-2012",        "hypoxiaSignaling-2006","stamlabTFs-2012",
                "recon202",             "APID",                 "BioGrid",
                "bhf-ucl",              "ChEMBL",               "DIP",
                "HPIDb",                "InnateDB",             "IntAct",
                "mentha",               "MPIDB",                "MatrixDB",
                "MINT",                 "Reactome",             "Reactome-FIs",
                "STRING",               "BIND",                 "Interoporc",
                "I2D-IMEx",             "InnateDB-IMEx",        "MolCon",
                "UniProt",              "MBInfo",               "BindingDB",
                "VirHostNet",           "Spike",                "BAR")

#----------------------------------------------------------------------------------------------------
javascript <- function()
{
    text <-  '<script>
var geneTextBox  = $("#genes");
var providerSelector = $("#providers");
var goButton  = $("#goButton");

assessGoButtonReadiness = function() {
   var gotPossibleGeneName = geneTextBox.val().length > 1;
   var gotProvider = providerSelector.val() != null
   if(gotPossibleGeneName & gotProvider)
      goButton.prop("disabled", false);
   else
      goButton.prop("disabled", true);
   } // assessGoButtonReadiness

geneInputBoxReader = function(e) {
   console.log("e.keyCode: " + e.keyCode);
   console.log("length: " + geneTextBox.val().length);
   assessGoButtonReadiness();
   }

providersReader = function() {
   assessGoButtonReadiness();
   }

 $(document).ready(function() {
    geneTextBox  = $("#genes")
    providerSelector = $("#providers")
    goButton  = $("#goButton");
    goButton.prop("disabled", true);

    geneTextBox.keydown(geneInputBoxReader);
    providerSelector.change(providersReader);

    $("#table").on("click", "tr", function () {
        console.log("table clicked!");
        //if($(this).hasClass("rowsSelected"))
        $("tr.rowsSelected").removeClass("rowsSelected")
        $(this).addClass("rowsSelected");
        var Aname = $("td", this).eq(0).text();
        var Bname = $("td", this).eq(1).text();
        var type  = $("td", this).eq(2).text();
        var pmid =  $("td", this).eq(3).text();
        var A_id = $("td", this).eq(6).text();
        var B_id = $("td", this).eq(7).text();
        //saveButton = $("#saveInteractionButton");
        //console.log("saveButton: " + saveButton);
        //saveButton.click(function() console.log("save!"));
        console.log("A_id: " + A_id);
        window.pmid = pmid;
        console.log("pmid ? " + pmid)

        $("#hiddenPmidDiv").html(pmid);
        $("#hiddenPmidDiv").hide();
        Shiny.onInputChange("hiddenPmidDiv", pmid);

        $("#hiddenGeneADiv").html(Aname);
        $("#hiddenGeneADiv").hide();
        Shiny.onInputChange("hiddenGeneADiv", Aname);

        $("#hiddenGeneBDiv").html(Bname);
        $("#hiddenGeneBDiv").hide();
        Shiny.onInputChange("hiddenGeneBDiv", Bname);

        $("#hiddenGeneA_id_Div").html(A_id);
        $("#hiddenGeneA_id_Div").hide();
        Shiny.onInputChange("hiddenGeneA_id_Div", A_id);

        $("#hiddenGeneB_id_Div").html(B_id);
        $("#hiddenGeneB_id_Div").hide();
        Shiny.onInputChange("hiddenGeneB_id_Div", B_id);
        })});
</script>'

    return(text)

} # javascript
#----------------------------------------------------------------------------------------------------
rowSelectableDataTable <- function(outputId, ...)
{
     origStyle<- c( 
         '<script src="shared/datatables/js/jquery.dataTables.min.js"></script>',
         '<script class="shiny-html-output" 
                  src= "/js/DTbinding.js"></script>',
         '<link rel = "stylesheet", 
                type = "text/css", 
                href = "shared/datatables/css/DT_bootstrap.css"></link>',
         '<style type="text/css">
                .rowsSelected td{
                background-color: rgb(204,255,204) !important}  </style>',
         '<style type="text/css"> .selectable div table tbody tr{
                cursor: hand; cursor: pointer;}</style>',
         '<style type="text/css"> .selectable div table tbody tr td{
                -webkit-touch-callout: none;
                -webkit-user-select: none;
                -khtml-user-select: none;
                -moz-user-select: none;
                -ms-user-select: none;
                user-select: none;} </style>',
         '<style type="text/css">
                #myTable tfoot {display:table-header-group;}
               th, td { white-space: nowrap; }
         </style>')
     
     tagList(
         singleton(
             tags$head(HTML(origStyle))
             ),
         div(id = outputId, class = "shiny-datatable-output selectable")
         )

} # rowSelectableDataTable
#----------------------------------------------------------------------------------------------------
shinyUI <- fluidPage(
   HTML(javascript()),
   headerPanel("RefNet (Homo sapiens)"),
       sidebarPanel(width=3,
          selectizeInput(inputId='providers', label='providers',
                         choices = providers, multiple = TRUE,
                         selected = NULL), # "mentha"),
                         #selected = initial.providers),
                         #selected=GLOBALS[["initial.providers"]]),
          textInput("genes", "Genes:", value=""),
          actionButton("goButton", "Find interactions"),
          #actionButton("saveInteractionButton", "Save Interaction"),
          htmlOutput("hiddenPmidDiv"),
          htmlOutput("hiddenGeneADiv"),
          htmlOutput("hiddenGeneBDiv"),
          htmlOutput("hiddenGeneA_id_Div"),
          htmlOutput("hiddenGeneB_id_Div")
          ),
       mainPanel(
          tabsetPanel(
               tabPanel('interactions', rowSelectableDataTable(outputId="table")),
               tabPanel("Pubmed Abstract", htmlOutput("pubmedAbstract")),
               tabPanel("A", htmlOutput("geneA")),
               tabPanel("B", htmlOutput("geneB"))
               #tabPanel("Saved Interactions", rowSelectableDataTable(outputId="selectedInteractionsTable"))
             )
          )
       ) # uiWidgets

#----------------------------------------------------------------------------------------------------

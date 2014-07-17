setGeneric("explore", signature="object",
             function(object,
                      id=NA,
                      species=NA,
                      speciesExclusive=TRUE,
                      type=NA,
                      provider=NA,
                      detectionMethod=NA,
                      publicationID=NA,
                      quiet=TRUE,
                      port=9999)
           standardGeneric ("explore"))
#-------------------------------------------------------------------------------
setMethod("explore", "RefNet",

   function (object, 
             id=NA,
             species=NA,
             speciesExclusive=TRUE,
             type=NA,
             provider=NA,
             detectionMethod=NA,
             publicationID=NA,
             quiet=TRUE,
             port=9999){
       printf("starting shiny...");
       appDir <- system.file(package="RefNet", "apps")
       stopifnot(file.exists(appDir))
       stopifnot(file.exists(file.path(appDir, "server.R")))
       stopifnot(file.exists(file.path(appDir, "ui.R")))
       shiny::runApp(appDir=appDir, port=port)
       })

#-------------------------------------------------------------------------------


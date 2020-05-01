#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
library(ggalluvial)

source("enrichment.R")

lib.list <- dir(path = "../../Enrichment Libs", "*")

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Enrichment Report"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            selectInput("dt", "Data Type:", c("SOMA", "Metabolon", "RNA-Seq")),
            fileInput("exp_names", "Expressed Features", multiple = F),
            fileInput("obs_names", "Selected Features", multiple = F),
            selectInput("lib1", "Library:", lib.list),
            sliderInput("k", "k:", min=1, max=20, value=5),
            sliderInput("rank", "Rank:", min=1, max=500, value=1),
            downloadButton("download_res", label = "Download as Table")
        ),

        # Show a plot of the generated distribution
        mainPanel(
            textOutput("exp.annotated"),
            textOutput("obs.annotated"),
            tabsetPanel(
                tabPanel("Plot", plotOutput("barplot")),
                tabPanel("Table", tableOutput("result.tab")),
                tabPanel("Annotation", plotOutput("alluvial"))
            )
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    expressedLabel <- reactive({
        if (!is.null(input$exp_names)) {
            exp.names <- read.delim(input$exp_names$datapath, header=F)
            exp.names$V1
        } else {
           NULL
        }
    })
    
    readObserved <- reactive({
        if (!is.null(input$obs_names)) {
            obs <- read.delim(input$obs_names$datapath, sep=",")
            obs
        } else {
            NULL
        }
    })
    
    observedLabel <- reactive({
        obs.names <- readObserved()
        if (!is.null(obs.names)) {
            obs.names <- as.vector(obs.names[obs.names$ranks <= input$rank, "Column"])
            if(input$dt == "RNA-Seq"){
                obs.names <- gene.ids2symbols(obs.names)
                obs.names$hgnc_symbol
            } else {
                obs.names
            }
        } else {
            NULL
        }
    })
    
    expressedAnnotate <- reactive({
        exp.names <- expressedLabel()
        if (!is.null(exp.names)){
            if(input$dt == "RNA-Seq"){
                exp.names <- gene.ids2symbols(exp.names)$hgnc_symbol
                # exp <- gene.annotate.from(exp.names, paste0("../../Enrichment Libs/", input$lib1))
                exp <- PANTHER.annotate(exp.names)
            } else if(input$dt == "Metabolon"){
                exp <- met.annotate.from(exp.names, paste0("../../Enrichment Libs/", input$lib1))
            }
            exp
        } else {
            NULL
        }
        
    })
    
    observedAnnotate <- reactive({
        obs.names <- observedLabel()
        exp.names <- expressedAnnotate()
        if (!is.null(obs.names)){
            if(input$dt == "RNA-Seq"){
                # obs <- gene.annotate.from(obs.names, paste0("../../Enrichment Libs/", input$lib1))
                obs <- PANTHER.annotate(obs.names)
            } else if(input$dt == "Metabolon"){
                obs <- met.annotate.from(obs.names, paste0("../../Enrichment Libs/", input$lib1))
            } else if (input$dt == "SOMA"){
                obs <- soma.annotate.from(obs.names, paste0("../../Enrichment Libs/", input$lib1))
            }
            obs
        }
        else {
            NULL
        }
    })

    datasetEnrich <- reactive({
        if (!is.null(input$exp_names) & !is.null(input$obs_names)){
            exp <- expressedAnnotate()$PATHWAY
            obs <- observedAnnotate()$PATHWAY
            results <- enrich.fishers(exp[exp != "Unannotated"], obs[obs != "Unannotated"])
            results
        } else {
            NULL
        }
    })
    
    datasetAsNumeric <- reactive({
        res <- datasetEnrich()
        res$exp.count <- res$exp.freq
        res$obs.count <- res$obs.freq
        if (!is.null(res)){
            obs <- readObserved()
            frac <- length(obs[obs$ranks <= input$rank, "Column"]) / sum(res$exp.count) #length(obs$Column)
            res$exp.freq <- res$exp.count * frac
            res$obs.freq <- res$obs.count / res$exp.freq
            colnames(res) <- c("Pathway", "Expected", "Observed-Expected Ratio", "p-value (Fisher's)", 
                                "Adjusted p-value (FDR)", "Pathway Size",
                                "Observed")
        }
        res[,c(1, 6, 7, 2, 3, 4, 5)]
    })
    
    output$barplot <- renderPlot({
        exp.counts = rep(0, input$k)
        obs.counts = rep(0, input$k)
        p.vals = rep(0, input$k)
        paths <- rep("N/A", input$k)
        
        results <- datasetAsNumeric()
        
        if (!is.null(results)){
            exp.counts <- results$`Pathway Size`[1:input$k]
            obs.counts <- results$`Observed`[1:input$k]
            p.vals <- results$`p-value (Fisher's)`[1:input$k]
            paths <- results$Pathway[1:input$k]
        }
        
        plottables <- data.frame(pathw = rep(paths, 2), 
                                 counts=c(exp.counts, obs.counts), 
                                 observed=c(rep(0, length(exp.counts)), rep(1, length(obs.counts))),
                                 p=c(p.vals, p.vals))
        
        bp <- ggplot(data=plottables, aes(x=pathw, y=counts, fill=as.factor(observed))) +
                geom_bar(stat="identity") +
                geom_text(aes(label=round(p, digits = 4)), vjust=1.6, color="white", size=3.5) +
                theme(axis.text.x = element_text(angle = 60, hjust = 1))
                # position=position_dodge()
                #theme_minimal() +
                #coord_flip() 
        bp
        #plottables <- matrix(data=c(exp.counts[1:input$k], obs.counts[1:input$k]), nrow = 2)
        #bp = barplot(plottables, xlab = "Pathway", col=c("red", "grey"), names.arg = paths[1:input$k], las=2, beside = T)
        #text(bp, plottables[1,]+0.1, labels = round(p.vals, digits = 2))
    })
    
    output$result.tab <- renderTable(
        results <- datasetAsNumeric()
    )
    
    output$alluvial <- renderPlot({
        exp.ids <- expressedLabel()
        exp.syms <- gene.ids2symbols(exp.ids)$hgnc_symbol
        paths <- expressedAnnotate()$PATHWAY
        
        f1 <- length(exp.syms[!is.na(exp.syms)])
        f2 <- length(paths[paths != "Unannotated"])
        f3 <- f1 - f2
        f4 <- length(exp.syms[is.na(exp.syms)])
        f5 <- length(exp.ids) - length(exp.syms)
        
        annots <- data.frame(Ids="Valid", Symbols="Annotated", Pathways="Annotated", Freq=f2)
        annots <- rbind(annots, data.frame(Ids="Valid", Symbols="Annotated", Pathways="Unannotated", Freq=f3))
        annots <- rbind(annots, data.frame(Ids="Valid", Symbols="Unannotated", Pathways="Annotated", Freq=0))
        annots <- rbind(annots, data.frame(Ids="Valid", Symbols="Unannotated", Pathways="Unannotated", Freq=f4))
        annots <- rbind(annots, data.frame(Ids="Invalid", Symbols="Annotated", Pathways="Annotated", Freq=0))
        annots <- rbind(annots, data.frame(Ids="Invalid", Symbols="Annotated", Pathways="Unannotated", Freq=0))
        annots <- rbind(annots, data.frame(Ids="Invalid", Symbols="Unannotated", Pathways="Annotated", Freq=0))
        annots <- rbind(annots, data.frame(Ids="Invalid", Symbols="Unannotated", Pathways="Unannotated", Freq=f5))
        
        ap <- ggplot(data = annots, aes(y = Freq,
                   axis1 = Ids, axis2 = Symbols, axis3 = Pathways)) +
            geom_alluvium(aes(fill = Pathways),
                          width = 0, knot.pos = 0, reverse = FALSE) +
            guides(fill = FALSE) +
            geom_stratum(width = 1/8, reverse = FALSE) +
            geom_text(stat = "stratum", infer.label = TRUE, reverse = FALSE) +
            scale_x_continuous(breaks = 1:3, labels = c("Ids", "Symbols", "Pathways")) +
            coord_flip() +
            ggtitle("Annotation Breakdown")
        ap
    })
    
    # Downloadable csv of selected dataset ----
    output$download_res <- downloadHandler(
        filename = function() {
            exp <- expressedAnnotate()$PATHWAY
            annotated <- 0
            if(!is.null(exp)){
                annotated = sum(exp != "Unannotated")
            }
            paste(input$lib1, "_enrich_", annotated, ".csv", sep = "")
        },
        content = function(file) {
            write.csv(datasetAsNumeric(), file, row.names = FALSE)
        }
    )
    
    output$exp.annotated <- renderText({
        exp <- expressedAnnotate()$PATHWAY
        annotated <- 0
        exp.counts <- 0
        if(!is.null(exp)){
            annotated = sum(exp != "Unannotated")
            exp.counts = length(as.data.frame(table(exp))$exp)
        }
        totals <- expressedLabel()
        sprintf("Expressed:\n\tAnnotated: %i/%i\n\tUnique Pathways: %i", annotated, length(totals), exp.counts)
    })
    
    output$obs.annotated <- renderText({
        obs <- observedAnnotate()$PATHWAY
        annotated <- 0
        obs.counts <- 0
        if(!is.null(obs)){
            annotated = sum(obs != "Unannotated")
            obs.counts = length(as.data.frame(table(obs))$obs)
        }
        totals <- observedLabel()
        sprintf("Observed:\n\tAnnotated: %i/%i\n\tUnique Pathways: %i", annotated, length(totals), obs.counts)
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

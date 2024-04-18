#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# Load library
library(shiny)
library(ggplot2)
library(plotly)
library(dplyr)
library(tidyr) 
library(corrplot)
library(reshape2)
library(heatmaply)
library(stringr)

# load data
icr_data <- read.csv("/Users/vivek/Desktop/ICR_Final/RshinyApp/data/cICR1495_MatPatOg.csv", stringsAsFactors = FALSE,header=F,sep=",")
colnames(icr_data) <- c("ICR_name","Chromosome:start_position-end_position","maternal/paternal")
genes_data <- read.delim("/Users/vivek/Desktop/ICR_Final/RshinyApp/data/AllGenes.txt", stringsAsFactors = FALSE)
xpr_data <- read.csv("/Users/vivek/Desktop/ICR_Final/RshinyApp/data/xpr_lcpm.tsv", sep = ",", header = TRUE, row.names = 1)
exon_data <- read.delim("/Users/vivek/Desktop/ICR_Final/RshinyApp/data/AllGenesExons.tsv", header = FALSE, col.names = c("Region", "Exon_Start", "Exon_End"))
#asm_data_1 <- read.csv("/Users/vivek/Desktop/ICR_Final/RshinyApp/data/21424B.asm", header = TRUE,sep="\t") %>% mutate(Source = "21424B.asm")
#asm_data_2 <- read.csv("/Users/vivek/Desktop/ICR_Final/RshinyApp/data/21423L.asm", header = TRUE,sep="\t") %>% mutate(Source = "21423L.asm")
#asm_data_3 <- read.csv("/Users/vivek/Desktop/ICR_Final/RshinyApp/data/21424L.asm", header = TRUE,sep="\t") %>% mutate(Source = "21424L.asm")
#asm_data_4 <- read.csv("/Users/vivek/Desktop/ICR_Final/RshinyApp/data/21562L.asm", header = TRUE,sep="\t") %>% mutate(Source = "21562L.asm")
#asm_data_5 <- read.csv("/Users/vivek/Desktop/ICR_Final/RshinyApp/data/22288B.asm", header = TRUE,sep="\t") %>% mutate(Source = "22288B.asm")
#asm_data_6 <- read.csv("/Users/vivek/Desktop/ICR_Final/RshinyApp/data/22316B.asm", header = TRUE,sep="\t") %>% mutate(Source = "22316B.asm")

# Pre-process data
genes_data$hg38.knownGene.chrom <- sub("_.*", "", genes_data$hg38.knownGene.chrom)
genes_data$hg38.knownGene.chrom <- sub("^chr", "Chr", genes_data$hg38.knownGene.chrom)
#combined_asm_data <- bind_rows(asm_data_1, asm_data_2, asm_data_3, asm_data_4, asm_data_5, asm_data_6)
#combined_asm_data$Chr <- sub("_.*", "",combined_asm_data$Chr)

# Split the chromosome:start_position-end_position into separate columns for easier processing.
icr_data <- icr_data %>%
    separate(`Chromosome:start_position-end_position`, into = c("Chromosome", "Positions"), sep = ":") %>%
    separate(Positions, into = c("Start", "End"), sep = "-") %>%
    mutate(Start = as.numeric(Start), End = as.numeric(End))

icr_data$dropdown_display <- paste(icr_data$ICR_name, "(", icr_data$Chromosome, ":", icr_data$Start, "-", icr_data$End, ")", sep = "")


# Define UI for application that draws a histogram
ui <- fluidPage(
    titlePanel("Imprinted Genes Viewer"),
    sidebarLayout(
        sidebarPanel(
            selectInput("selected_icr", "Select ICR:", choices = icr_data$ICR_name),
            sliderInput("range", "Range around ICR:", min = 250, max = 1000, value = 1000, step = 250, post = "Kb"),
            textOutput("gene_count"),
            textOutput("mRNA_seq_gene_count"),
            selectInput("gene_list", "Select Gene:", choices = NULL),
            width = 2
        ),
        mainPanel(
            plotlyOutput("chromosome_plot",height = "800px", width = "100%"),
            #selectInput("gene_list", "Select Gene:", choices = NULL),  # Add a select input for gene names
            plotlyOutput("correlation_plot",height = "800px", width = "100%"),
            dataTableOutput("gene_expression_table"),
            dataTableOutput("snp_table"),  
            uiOutput("gene_image") 
        )
    )
)

server <- function(input, output, session) {
    selected_icr_data <- reactive({
        selected_data <- icr_data[icr_data$ICR_name == input$selected_icr, ]
        if (nrow(selected_data) == 0) {
            return(data.frame(ICR_name = NA, Chromosome = NA, Start = NA, End = NA, stringsAsFactors = FALSE))
        }
        selected_data
    })
    
    filtered_genes_data <- reactive({
        selected_icr <- icr_data[icr_data$ICR_name == input$selected_icr, ]
        window_start <- selected_icr$Start - (input$range * 1000)
        window_end <- selected_icr$End + (input$range * 1000)
        filtered_genes <- genes_data %>%
            filter(hg38.knownGene.chrom == selected_icr$Chromosome,
                   hg38.knownGene.txStart <= window_end,
                   hg38.knownGene.txEnd >= window_start) %>%
            distinct(hg38.kgXref.geneSymbol, hg38.knownGene.chrom) %>%
            select(hg38.kgXref.geneSymbol, hg38.knownGene.chrom) %>%
            filter(hg38.kgXref.geneSymbol %in% rownames(xpr_data))
        filtered_genes
    })
    
    gene_expression_count <- reactive({
        nrow(filtered_genes_data())
    })
    
    output$mRNA_seq_gene_count <- renderText({
        count <- gene_expression_count()
        paste("Total number of genes with mRNA-seq expression:", count)
    })
    
    output$chromosome_plot <- renderPlotly({
        icr_data$selected <- icr_data$ICR_name == input$selected_icr
        icr_data$Chromosome <- factor(icr_data$Chromosome, levels = c(paste0("Chr", 1:22), "ChrX", "ChrY", "ChrMT"))
        p <- ggplot(icr_data, aes(x = Start, y = Chromosome, text = ICR_name, customData=ICR_name)) +
            geom_point(aes(color = selected), size = 2) +
            scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
            labs(title = "Chromosome Plot", x = "Position", y = "Chromosome") +
            theme_minimal()
        ggplotly(p, tooltip = "text", source="chromosomePlot") %>% 
            layout(hovermode = "closest", hoverlabel = list(bgcolor = "white", font = list(size = 12), bordercolor = "black", align = "left"))
    })
    
    observeEvent(event_data("plotly_click", source = "chromosomePlot"), {
        click_data <- event_data("plotly_click", source = "chromosomePlot")
        if (!is.null(click_data)) {
            point_number <- click_data[["pointNumber"]]
            if (!is.null(point_number) && point_number + 1 <= nrow(icr_data)) {
                clicked_icr_name <- icr_data$ICR_name[point_number + 1]  # Adjust index if needed
                if (clicked_icr_name %in% icr_data$ICR_name) {
                    updateSelectInput(session, "selected_icr", selected = clicked_icr_name)
                }
            }
        }
    })
    
    output$gene_table <- renderDataTable({
        filtered_genes_data()
    })
    
    output$gene_count <- renderText({
        count <- nrow(filtered_genes_data())
        paste("Total number of genes within", input$range, "Kb upstream and downstream:", count)
    })

   
    output$correlation_plot <- renderPlotly({
        filtered_genes <- filtered_genes_data()
        selected_genes <- filtered_genes$hg38.kgXref.geneSymbol
        
        if (length(selected_genes) > 1) {
            gene_data <- xpr_data[selected_genes, , drop = FALSE]
            correlation_matrix <- abs(cor(t(gene_data)))
            col_names <- colnames(gene_data)
            
            # Use heatmaply to create the heatmap with dendrograms
            heatmaply(
                gene_data,
                colors = colorRampPalette(c( "white", "#542788"))(256),
                show_dendrogram = TRUE,
                dendrogram_height = 0.2
            )
        } else {
            plot.new()
            text(0.5, 0.5, "Not enough genes for correlation plot", cex = 1.5)
        }
    })
    
    observe({
        gene_names <- filtered_genes_data()$hg38.kgXref.geneSymbol
        updateSelectInput(session, "gene_list", choices = gene_names)
    })

    output$gene_expression_table <- renderDataTable({
        selected_gene <- input$gene_list
        if (!is.null(selected_gene) && selected_gene %in% rownames(xpr_data)) {
            gene_expression <- t(data.frame(xpr_data[selected_gene, , drop = FALSE]))
            colnames(gene_expression) <- c("xpr_lcpm")
            gene_expression
        } else {
            data.frame()  # Return an empty data frame if no gene is selected or the gene is not in xpr_data
        }
    })
    
    
}


# Run the application 
shinyApp(ui = ui, server = server)

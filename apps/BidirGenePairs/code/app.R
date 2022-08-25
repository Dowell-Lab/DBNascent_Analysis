library(shiny) ##shiny interactions
library(ggplot2) ##plotting package
library(cowplot) ##plots + arrange plots etc
library(RColorBrewer) ## color palettes
library(dplyr) ## for the R pipes
library(tidyr) ## for tidying the data.frames
library(data.table) ## load files into R faster
library(scales) ## representing values with varying scales
library(DT) ##interactive datatables from a java script library 

#Import data
tpms <- data.table::fread("~/Desktop/cotranscription_analysis_paper/BidirGenePairs/data/chr21_genes_bidirectionals_tpm.tsv.gz")
metadata <- read.csv("~/Desktop/cotranscription_analysis_paper/BidirGenePairs/data/20211026_db_output.csv")
tissues <- read.csv("~/Desktop/cotranscription_analysis_paper/BidirGenePairs/data/human_cell_types.csv")
chr21_corr_pairs <- data.table::fread("~/Desktop/cotranscription_analysis_paper/BidirGenePairs/data/chr21_corrs_keep_vs_remove0s_pairs2mb.tsv.gz")

##create master metadata-tissue table 
metadata_celltype <- merge(metadata, tissues, by="cell.type", all=TRUE)

## define number of samples
nsamples <- 880

## select colors to use
# Define the number of colors you want
nb.cols <- 20
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)

##get transcript names
gene_names <- unique(chr21_corr_pairs$transcript_1.x)
bidir_names <- unique(chr21_corr_pairs$transcript_2.x)

pair_ids <- unique(chr21_corr_pairs$pair_id)

#Functions
plot_coefficients <- function(corrs_compared, gene_id, bidir_id){
  
  r2mb <- ggplot(corrs_compared, 
                 aes(coefficient.x, coefficient.y)) +
    geom_point(shape=21, 
               size=2, 
               fill='maroon2',
               color='gray40', 
               alpha=0.05) +
    ggtitle("Chromosome 21 Pairs \n within 2MB Distance") +
    xlab("Pearson's R \n Keep 0s")+
    ylab("Pearson's R \n Remove samples with 0s") + 
    theme_cowplot() + 
    scale_x_continuous(limits =c(-1,1)) +
    scale_y_continuous(limits =c(-1,1)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.title = element_text(hjust = 0.5),
          title = element_text(size = 16), 
          axis.title = element_text(size = 14), 
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12)) +
    geom_point(data = subset(corrs_compared,
                             transcript_1.x == gene_id &
                               transcript_2.x == bidir_id),
               aes(coefficient.x, coefficient.y),
               color='black',
               fill='maroon4',
               size=4,
               shape=21) +
    geom_hline(yintercept=0, 
               linetype="dashed", 
               color = "gray80") +
    geom_vline(xintercept=0, 
               linetype="dashed",
               color = "gray80")
  
  print(r2mb)
  
}

filter_transform_tpms <- function(normalized_counts, 
                                  nsample, gene_id, 
                                  bidir_id, metadata_df){
  #get sample numbers
  nsamples <- nsample
  
  #convert counts DT to data.frame
  normalized_counts_df <- as.data.frame(normalized_counts)
  
  #filter gene and bidirectional transcripts
  gene_filtered <- subset(normalized_counts_df, GeneID == as.character(gene_id))
  bidir_filtered <- normalized_counts_df[normalized_counts_df$GeneID %in% bidir_id,]
  
  #combine
  filtered_pair <- rbind(gene_filtered, bidir_filtered)
  rownames(filtered_pair) <- filtered_pair$GeneID
  
  #convert the count to long form
  filtered_pair_df <- as.data.frame(log(t(filtered_pair[c(6:(nsample+5))])+1, base=10))#as.data.frame(t(filtered_pair[c(6:(nsample+5))]))
  colnames(filtered_pair_df) <- c('transcript_1','transcript_2')
  filtered_pair_df$srz <- rownames(filtered_pair_df)
  
  #merge with metadata
  filtered_pair_metadata <- merge(filtered_pair_df, 
                                  metadata_df, 
                                  by.x = 'srz',
                                  by.y = 'srz')
  
  #Log transform counts and convert 0s to NAs for correlation plotting
  filtered_pair_log <- as.data.frame(log(t(filtered_pair[c(6:(nsample+5))])+1, base=10))
  colnames(filtered_pair_log) <- c('transcript_1','transcript_2')
  filtered_pair_log$srz <- rownames(filtered_pair_log)
  filtered_pair_log_NAs <- dplyr::na_if(filtered_pair_log, 0)
  
  #Merge the long counts with metadata
  filtered_pair_log_metadata <- merge(filtered_pair_log_NAs,
                                      metadata_df, 
                                      by.x = 'srz',
                                      by.y = 'srz')
  
  #return both TPMs and log transformed TPMs with metadata
  return(list(filtered_pair_metadata, filtered_pair_log_metadata))
  
}

plot_correlations <- function(transformed_tpms, datatable_pairs,
                              gene_id, bidir_id, cols,
                              remove0s=TRUE){
  
  #tpms_keep_0s <- transformed_tpms[[1]]
  #tpms_remove_0s <- transformed_tpms[[2]]
  
  #subset the table with both 0 and NA correlations
  corr_summary <- subset(datatable_pairs, 
                         transcript_1.x == gene_id &
                           transcript_2.x == bidir_id)
  
  #get transcript ids
  transcript_1id <- corr_summary$transcript_2.x 
  transcript_2id <- corr_summary$transcript_1.x 
  
  if(remove0s==TRUE){
    
    #get coefficients and p-values
    coeff <- corr_summary$coefficient.y
    adj <- corr_summary$adj_p_BH.y
    
  } else {
    coeff <- corr_summary$coefficient.x
    adj <- corr_summary$adj_p_BH.x
  }
  
  #get coefficients and p-values
  #coeff_x <- corr_summary$coefficient.x
  #adj_x <- corr_summary$adj_p_BH.x
  #coeff_y <- corr_summary$coefficient.y
  #adj_y <- corr_summary$adj_p_BH.y
  
  #scatter plot with line fit
  corr_plot <- ggplot(transformed_tpms, 
               aes(x=transcript_2,
                   y=transcript_1,
                   fill=tissue)) + 
    geom_point(size=3,shape=21, color='gray80')+
    theme_cowplot() + 
    geom_smooth(method='lm',
                formula=y~x,
                color='gray40',
                fill = 'gray60') +
    labs(title = paste("R = ",
                       round(as.numeric(coeff), 3),
                       ",", 
                       "Adjusted P = ",
                       round(as.numeric(adj), 3)),
         x = bquote(.(transcript_1id) ~ log[10] ~ TPM), 
         y = bquote(.(transcript_2id) ~ log[10] ~ TPM)) + 
    scale_fill_manual(name = "Tissues", values = cols) +
    theme(plot.title = element_text(hjust = 0.5,
                                    size = 16),
          axis.title = element_text(size = 14), 
          axis.text = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10),
          axis.title.x = element_text(colour = "#5e3c99"),
          axis.title.y = element_text(colour = "#e66101"),
          panel.background = element_rect(fill = "transparent"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          legend.key = element_rect(colour = "transparent", 
                                    fill = "transparent"),
          legend.background = element_rect(fill = "transparent"), 
          legend.box.background = element_rect(fill = "transparent",
                                               colour = "transparent"))   
  
  print(corr_plot)
  
}

plot_boxplots <- function(transformed_tpms, datatable_pairs, 
                          gene_id, bidir_id,cols,
                          transcript_id='transcript_1'){
  
  #subset the table with both 0 and NA correlations
  corr_summary <- subset(datatable_pairs, 
                         transcript_1.x == gene_id &
                           transcript_2.x == bidir_id)
  
  #get transcript ids
  transcript_1id <- corr_summary$transcript_2.x 
  transcript_2id <- corr_summary$transcript_1.x 
  
  #total number of samples with transcripts
  if(transcript_id =='transcript_1'){
    ##transcript 1 summary
    summ_transcript <- transformed_tpms[c('transcript_1',
                                           'tissue')] %>% 
      filter(!is.na(transcript_1)) %>%
      group_by(tissue) %>% 
      summarise(n = n(),
                transcript_1 = max(transcript_1))
    
    transcript_name <- corr_summary$transcript_1.x 
    
    id_col <- "#e66101"

  } else {
    ##transcript 2 summary
    summ_transcript <- transformed_tpms[c('transcript_2',
                                          'tissue')] %>% 
      filter(!is.na(transcript_2)) %>%
      group_by(tissue) %>% 
      summarise(n = n(),
                transcript_2 = max(transcript_2))
    
    transcript_name <- corr_summary$transcript_2.x 
    
    id_col <- "#5e3c99"
  }
  
  #plot box plot for transcript
  transcript_boxplot <- ggplot(transformed_tpms,
               aes(x=tissue, 
                   y=get(transcript_id), 
                   fill=tissue)) +
    geom_boxplot(color='gray80') + 
    geom_text(aes(label = n), size=4, data = summ_transcript) +
    ggtitle(" ") +
    xlab(" ") +
    ylab(bquote(.(transcript_name) ~ log[10] ~ TPM)) +
    scale_fill_manual(name = "Tissues",
                      values = cols, drop = FALSE) +
    theme_cowplot() +
    theme(plot.title = element_text(hjust = 0.5,
                                    size = 16),
          axis.title = element_text(size = 12), 
          axis.text = element_text(size = 14),
          axis.text.x = element_text(size = 14, 
                                     angle = 45,
                                     hjust=0.5, 
                                     vjust=0.5),
          axis.title.y = element_text(colour = id_col),
          panel.background = element_rect(fill = "transparent"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          legend.position = "none")
  
  print(transcript_boxplot)

}

# UI code
ui <- fluidPage(
  fluidRow(
    column(12, dataTableOutput("correlations_table"))
  ),
  fluidRow(
    column(6, selectInput("gene", "Gene ID", choices = gene_names)),
    column(6, selectInput("bidir", "Bidirectional ID", choices = bidir_names))
  ),
  fluidRow(
    column(6, plotOutput("coefficients")),
    column(6, tableOutput("correlations"))
  ),
  fluidRow(
    column(6, plotOutput("pair_with_zeros")),
    column(6, plotOutput("pair_without_zeros"))
    
  ),
  fluidRow(
    column(12, plotOutput("transcript1_tpms"))
  ),
  fluidRow(
    column(12, plotOutput("transcript2_tpms"))
  )
)

# Server code
server <- function(input, output, session) {
  output$correlations_table <- renderDataTable(chr21_corr_pairs[,c("pair_id",
                                                     "distance",
                                                     "coefficient.x",
                                                     "coefficient.y",
                                                     "r_diffrence")], 
                                               filter = 'top', 
                                               server = TRUE,
                                               options = list(
                                                       pageLength = 5))
  
  selected <- reactive(chr21_corr_pairs %>% filter(transcript_1.x == input$gene &
                                                   transcript_2.x == input$bidir))
  
  output$correlations <- renderTable(
    selected() %>% 
      select(pair_id, distance, r_diffrence) %>% gather(Name, Value)
  )
  
  output$coefficients <- renderPlot({
    plot_coefficients(chr21_corr_pairs, 
                      input$gene,
                      input$bidir)
  }, res = 96)
  
  output$pair_with_zeros <- renderPlot({
    
    filtered_tpms_list <- filter_transform_tpms(tpms,
                                                nsamples,
                                                input$gene,
                                                input$bidir,
                                                metadata_celltype)
    
    plot_correlations(filtered_tpms_list[[1]],
                      chr21_corr_pairs,
                      input$gene,
                      input$bidir,
                      mycolors,
                      remove0s = FALSE)
    
  }, res = 96)
  
  output$pair_without_zeros <- renderPlot({
    
    filtered_tpms_list <- filter_transform_tpms(tpms,
                                                nsamples,
                                                input$gene,
                                                input$bidir,
                                                metadata_celltype)
    
    plot_correlations(filtered_tpms_list[[2]],
                      chr21_corr_pairs,
                      input$gene,
                      input$bidir,
                      mycolors)
    
  }, res = 96)
  
  output$transcript1_tpms <- renderPlot({
    
    filtered_tpms_list <- filter_transform_tpms(tpms,
                                                nsamples,
                                                input$gene,
                                                input$bidir,
                                                metadata_celltype)
    
    plot_boxplots(filtered_tpms_list[[1]],
                      chr21_corr_pairs,
                      input$gene,
                      input$bidir,
                      mycolors)
    
  }, res = 96)
  
  output$transcript2_tpms <- renderPlot({
    
    filtered_tpms_list <- filter_transform_tpms(tpms,
                                                nsamples,
                                                input$gene,
                                                input$bidir,
                                                metadata_celltype)
    
    plot_boxplots(filtered_tpms_list[[2]],
                  chr21_corr_pairs,
                  input$gene,
                  input$bidir,
                  mycolors,
                  transcript_id = "transcript_2")
    
  }, res = 96)
  
}

shinyApp(ui, server)
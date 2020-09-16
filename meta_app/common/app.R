#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(GEOquery)
library(Biobase)
library(limma)
library(FoldGO)
library(MetaVolcanoR)
library(plotly)
library(dplyr)
library(sjmisc)
source("custom_draw.R")
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(factoextra))
source("metastudy_functions.R")
load("data/appdata.RData")


at_least <- function(num_vector, q_val){
  if(q_val < 0){
    return(num_vector <= q_val)
  }
  else if(q_val > 0){
    return(num_vector >= q_val)
  }
  else{
    return(NULL)
  }
}

#handling duplicate spots
median_dups <- function(logfc2, dupList){
  medians <- list()
  for(i in 1:length(dupList)){
    fc <- median(logfc2[as.vector(dupList[[i]])])
    medians[names(dupList)[i]] <- list(fc)
  }
  return(medians)
}
make_dup_list <- function(vnames){
  out <- list()
  i <- 1
  q_list <- vnames
  while(length(q_list > 1)){
    
    query <- vnames[i]
    
    q <- which(vnames == query)
    out[query] <- list(q)
    
    i <- i + 1
    
    q_list <- q_list[q]
    
  }
  #removing singles
  out <- out[which(as.vector(sapply(out, "length")) >= 2)]
  
  return(out)
  
}

put_na <- function(bool_set){
  if(bool_set){
    return(0)
  }
  else{
    return(NA)
  }
}
#requires a diffList to be loaded in the namespace and a list of gene names derived from query
pull_queried <- function(gene_ids){
  
  whole_set <- list()
  for(i in 1:length(diffList)){
    #picking apart valid and invalid ids
    valid_ids <- which(gene_ids %in% diffList[[i]]$Symbol)
    
    nvalid_ids <- which(!(gene_ids %in% diffList[[i]]$Symbol))
    
    bool_set <- gene_ids %in% diffList[[i]]$Symbol
    #adding NAs to unavailable ids
    log2fc <- sapply(bool_set, FUN = put_na)
    
    #dealing with duplicate entries by taking the median from the database
    whole_ids <- diffList[[i]]$Symbol[diffList[[i]]$Symbol %in% gene_ids]
    
    dups <- make_dup_list(whole_ids)
    
    if(length(dups) > 0){
      med_dups <- median_dups(logfc2 = diffList[[i]]$Log2FC, dupList = dups)
    
      #the duplicate location within the input id list
      dup_loc <- which(gene_ids %in% names(dups))
    
      log2fc[dup_loc] <- unlist(unname(med_dups))
    
      #removing valid_ids which are duplicates
      valid_ids <- valid_ids[which(!( valid_ids %in% dup_loc))]}
    
    
    log2fc[valid_ids] <-  diffList[[i]]$Log2FC[which(diffList[[i]]$Symbol %in% gene_ids[valid_ids])]
    names(log2fc) <- gene_ids
    
    whole_set[names(diffList)[[i]]] <- list(log2fc)
    
    
    
  }
  #browser()
  dfs <- lapply(whole_set, data.frame, stringsAsFactors = FALSE)
  binded.ids <- bind_cols(dfs)
  rownames(binded.ids) <- gene_ids
  colnames(binded.ids) <- names(diffList)
  return(binded.ids)
  
}


#Expects a list with named variables matching the ui tags
query <- function(varlist, smodel){
  out.ids <- list()
  print(varlist)
  model <- switch(smodel,
                  "Fisher's pvalue" = meta_comb@metaresult,
                  "Random Effects Modeling" = meta_degs_rem@metaresult)
  
  
  
  for(i in 1:length(varlist)){
    #getting the intersection of queries from the REM dataset
    if((varlist[i] != "") && smodel == "Random Effects Modeling"){
      name <- names(varlist)[i]
      out.var <- switch(name,"pval" = model["randomP"] <= as.numeric(varlist[i]),
                        "logfc2" = at_least(model["randomSummary"], as.numeric(varlist[i])),
                        "tRank" = model["rank"] <= as.numeric(varlist[i]),
                        "signcon" = at_least(model["signcon"], q_val = as.numeric(varlist[i]))
                        )
      
      out.ids[name] <- list(out.var)
    }
    
    
    
    
    else if((varlist[i] != "") && smodel == "Fisher's pvalue"){
      name <- names(varlist)[i]
      out.var <- switch(name, "pval" = model["metap"] <= as.numeric(varlist[i]),
                        "logfc2" = at_least(model["metafc"], as.numeric(varlist[i]))
                        )
      
      out.ids[i] <- list(out.var)
     
      
    }
    
  }
  dfs <- lapply(out.ids, data.frame, stringsAsFactors = FALSE)
  binded.ids <- bind_cols(dfs)
  rownames(binded.ids) <- model$Symbol
  #print(head(binded.ids))
  return(row.check(binded.ids))
}
row.check <- function(df){
  #rotate my data.frame bb
  df <- t(df)
  #output
  out <- c()
  #iterate and find the rows that are all true
  for(col in 1:length(df[1,])){
    if(all(df[,col])){
      out <- append(out, colnames(df)[col]) 
    }
  }
  return(out)
  
}



# Define UI for application 
ui <- fluidPage(
   
   # Application title
   titlePanel("Astroyeast MultiStress Explorer"),
   
   # Sidebar 
   sidebarLayout(
      sidebarPanel(
        helpText("helptext."),
         selectInput("var",
                     label = ("Choose a pvalue combination method"),
                     choices = c("Random Effects Modeling",
                                 "Fisher's pvalue"),
                     selected = "Random Effects Modeling"),
        
         textInput("gene", label = "Gene Level Cross Forest Plot", 
                   value = "Enter Common Gene Name"),
         helpText("choose threshold values to save genes for further analysis"),              
         fluidRow(
           column(5,
                  selectInput('smodel', 'Summary Model',c("","Random Effects Modeling","Fisher's pvalue"))),
           column(3,
                 textInput('pval',label = "P-value threshold")),
           column(3,
                 textInput('logfc2', label = "Summary Log2-Fold Change")),
           
           conditionalPanel(condition = "input.smodel == 'Random Effects Modeling'",
            column(3,
                    textInput('tRank', 'Number of Top Ranking Genes')),
           
           
            column(3, textInput("signcon", label = "sign Consistency"))
           )
         
          ),
         actionButton("save", label = "Save IDs"),
         helpText("Number of Queried Genes:"),
         verbatimTextOutput("table_summary")
         
        
        
       
      ),
        
      
      
      
      # Show a plot of the generated distribution
      
      
      mainPanel(
        
        tabsetPanel(type = "tabs",
          tabPanel("About", htmlOutput("about")),
          tabPanel("Summary Plot", 
            plotlyOutput("volcano"),
            plotOutput("forest")),
          
          tabPanel("Selected Gene Analysis", plotOutput("heat"), plotOutput("pca"))
          
          
          )
      )
        
        
        
      
   )
)

# Define server logic 
initial_txt <- function(string){
  return(string == "Enter Common Gene Name" || string == "")
}

server <- function(input, output) {
   
  #Rendering HTML Output for the ABOUT page 
  
  #Rendering main page graphs 
   output$volcano <- renderPlotly({
            switch(input$var, 
            "Random Effects Modeling" = ggplotly(meta_degs_rem@MetaVolcano),
            "Fisher's pvalue" = ggplotly(meta_comb@MetaVolcano))})
   
   

   
   output$forest <- renderPlot({if(!(initial_txt(input$gene))){
                                    
                                  draw_forest2(remres = meta_degs_rem,
                                            gene = input$gene,
                                            draw = "")}
                                            })
   
   #Saving and loading session data
   
   
   
   #query function based off of predefined parameters 
   reactive_query <- reactive({query(varlist = list(pval = input$pval, logfc2 = input$logfc2, tRank = input$tRank, signcon = input$signcon), input$smodel)})
   reactive_pull_query <- reactive({pull_queried(reactive_query())})
   
   output$table_summary <- renderText(length(reactive_query()))
   
   
  
   
   #When save ids is clicked, save them in session
   observeEvent(input$save, {
     reactive_data$ids <- reactive_query()
     reactive_data$logfc <- reactive_pull_query()
     reactive_data$pca <- reactive_pca()
     print(reactive_data$logfc)
     
   })
   
   #reactive dataset
   reactive_data <- reactiveValues(
     logfc = data.frame(),
     ids = vector(),
     pca = NULL
   )
   
   reactive_pca <- reactive({prcomp(na.omit(reactive_data$logfc))})
   
   
   
   output$heat <- renderPlot({Heatmap(as.matrix(na.omit(reactive_data$logfc)), heatmap_legend_param = list(title = "Log2FC"), show_row_names = FALSE)})
   output$pca <- renderPlot({fviz_pca_var(reactive_data$pca)})
   
   #print(reactive_pull_query())
   
}

# Run the application 
shinyApp(ui = ui, server = server)


library(shiny)
library(data.table)
library(shinydashboard)
library(shinymanager)
library(plyr)
library(e1071)
library(pROC)
library(dbscan)
library(DT)
library(linkcomm)
library(org.Hs.eg.db)
library(GSEABase)
library(GOstats)
source("function.R")
source("train.R")
source("clustering_LinkModule.R")
source("clustering_LinkPartition.R")
source("clustering_DBSCAN.R")
source("gotermenrichment.R")
source("LinkCluster.R")

inactivity <- "function idleTimer() {

var t = setTimeout(logout, 120000);
window.onmousemove = resetTimer; // catches mouse movements
window.onmousedown = resetTimer; // catches mouse movements
window.onclick = resetTimer;     // catches mouse clicks
window.onscroll = resetTimer;    // catches scrolling
window.onkeypress = resetTimer;  //catches keyboard actions

function logout() {
window.close();  //close the window
}
function resetTimer() {
clearTimeout(t);
t = setTimeout(logout, 120000);  // time is in milliseconds (1000 is 1 second)
}}
idleTimer();"

# data.frame with credentials info
credentials <- data.frame(
  user = c("1", "fanny", "victor", "benoit"),
  password = c("1", "azerty", "12345", "azerty"),
  # comment = c("alsace", "auvergne", "bretagne"), %>%
  stringsAsFactors = FALSE
)


# ui <- secure_app(head_auth = tags$script(inactivity),
ui <- dashboardPage(
      dashboardHeader(title = "Module identification system",titleWidth = 300),
      dashboardSidebar(
        sidebarMenu(
          menuItem("upload data",tabName = "a"),
          br(),
          menuItem("Clustering",tabName = "b"),
          br(),
          menuItem("Enrichment Analysis",tabName = "c"))
      ),
      dashboardBody(
        tabItems(
          tabItem(tabName = "a",
                  fluidRow(
                    box(
                      title = "input network data:",width = 6,solidHeader = TRUE,status = "primary",
                      fileInput("TF_miRNA","upload TF_miRNA network txt file",accept = ".txt"),
                      fileInput("TF_mRNA","upload TF_mRNA network txt file",accept = ".txt"),
                      fileInput("miRNA_mRNA","upload miRNA_mRNA network txt file",accept = ".txt"),
                      fileInput("mRNA_mRNA","upload mRNA_mRNA network txt file",accept = ".txt")
                    ),
                    box(
                      title = "input expression data:",width = 6,solidHeader = TRUE,status = "warning",
                      fileInput("miRNA","upload miRNA expression csv file",accept = ".csv"),
                      fileInput("TF","upload TF expression csv file",accept = ".csv"),
                      fileInput("mRNA","upload mRNA expression csv file",accept = ".csv")
                    )
                  ),
                  fluidRow(
                    box(
                      title = "input similarity data:",width = 6,solidHeader = TRUE,status = "danger",
                      fileInput("network_similarity","upload similarity_all csv file",accept = ".csv")
                      )
                    )
                  ),
          tabItem(tabName = "b",
                  fluidRow(
                    box(title = "Train network Represent Learning",width = 5,solidHeader = TRUE,status = "primary",
                        numericInput(inputId = "dim", label = "Training parameter dimension", value = 80, min = 10, max = 1000, step = 10, width = "60%"),
                        br(),
                        numericInput(inputId = "kmax", label = "Training parameter k-max", value = 50, min = 10, max = 1000, step = 10, width = "60%"),
                        br(),
                        actionButton(inputId = "train_button", label = "Submit! Go!", icon = icon("refresh"), style = "margin-right:20px"),height = 350
                        ),
                    box(title = "Train network Represent Learning output",width = 6,solidHeader = TRUE,status = "warning",
                        plotOutput("ROC",height = 285)
                        )
                    ),
                  fluidRow(
                    box(title = "Clustering",width = 5,solidHeader = TRUE,status = "danger",
                        "Choose the clustering method",
                        br(),
                        selectInput(inputId = "method", label = "  ", choices = c("LinkModule", "LinkPartition", "DBSCAN"), selected = FALSE),
                        "Click the button to perform clustering",
                        br(),br(),
                        actionButton(inputId = "cluster_button", label = "Submit Go!", icon = icon("refresh"), style = "margin-right:20px")
                    ),
                    box(title = "Clustering output",width = 6,solidHeader = TRUE,status = "info",
                        br(),
                        downloadButton(outputId = "download_result", label = "Download clustering result"),
                        br(),br(),
                        downloadButton(outputId = "download_target", label = "Download targetgene"),
                        br(),br()
                        # downloadButton(outputId = "download_target", label = "Download taget_gene")
                    )
                  )
                  ),
          tabItem(tabName = "c",
                  fluidRow(
                    box(title = "KEGG Enrichment Analysis",solidHeader = TRUE,status = "primary",width = 5,
                        # fileInput("targetgene","upload targetgene txt file",accept = ".txt"),
                        # "Click the button to perform go-terms analysis",br(),br(),
                        "upload target_gene txt file",br(),
                        fileInput("target"," ",accept = ".txt")
                        # actionButton(inputId = "goterm_button", label = "Submit! Go!", icon = icon("refresh"), style = "margin-right:20px")
                        ),
                    box(title = "Download KEGG pathways",solidHeader = TRUE,status = "danger",width = 5,"download KEGG pathways csv file",
                        br(),br(),
                        downloadButton(outputId = "download_kegg",label = "Download KEGG pathways"),
                        br(),br(),br()
                    )
                  ),
                    fluidRow(
                      box(title = "KEGG pathways Enrichment analysis",
                          dataTableOutput("summary"),width = 100)
                    )
                  ),
          tabItem(tabName = "d",tableOutput("net"))
          
    )
  )
)
# )

server <- function(input, output,session) {
  # result_auth <- secure_server(check_credentials = check_credentials(credentials))
  options(shiny.maxRequestSize=30*1024^2) 
  #####network data
  TF_miRNA <- reactive({
    inFile1 <- input$TF_miRNA
    if (is.null(inFile1)) return(NULL)
    read.csv(inFile1$datapath,header = FALSE,sep = "\t", stringsAsFactors = F)
  })
  TF_mRNA <- reactive({
    inFile2 <- input$TF_mRNA
    if (is.null(inFile2)) return(NULL)
    read.csv(inFile2$datapath,header = FALSE,sep = "\t", stringsAsFactors = F)
  })
  miRNA_mRNA <- reactive({
    inFile3 <- input$miRNA_mRNA
    if (is.null(inFile3)) return(NULL)
    read.csv(inFile3$datapath,header = FALSE,sep = "\t", stringsAsFactors = F)
  })
  mRNA_mRNA <- reactive({
    inFile4 <- input$mRNA_mRNA
    if (is.null(inFile4)) return(NULL)
    read.csv(inFile4$datapath,header = FALSE,sep = "\t", stringsAsFactors = F)
  })
  #######expression data
  miRNA <- reactive({
    inFile5 <- input$miRNA
    if (is.null(inFile5)) return(NULL)
    read.csv(inFile5$datapath,header = TRUE, row.names = 1, sep = ",",stringsAsFactors = F)
  })
  TF <- reactive({
    inFile6 <- input$TF
    if (is.null(inFile6)) return(NULL)
    read.csv(inFile6$datapath,header = TRUE, row.names = 1, sep = ",",stringsAsFactors = F)
  })
  mRNA <- reactive({
    inFile7 <- input$mRNA
    if (is.null(inFile7)) return(NULL)
    read.csv(inFile7$datapath,header = TRUE, row.names = 1, sep = ",",stringsAsFactors = F)
  })
  ###similarity data
  network_similarity <- reactive({
    inFile8 <- input$network_similarity
    if (is.null(inFile8)) return(NULL)
    read.csv(inFile8$datapath,header = T, sep = ",", row.names = 1)
  })
  
  target_gene <- reactive({
    inFile9 <- input$target
    if (is.null(inFile9)) return(NULL)
    read.table(inFile9$datapath,header = F, sep = "\t")
  })
  
  #####function######
  plot_fun_train <- eventReactive(input$train_button, {
    roc1 <- train(TF_miRNA(),TF_mRNA(),miRNA_mRNA(),mRNA_mRNA(),miRNA(),TF(),mRNA(),input$dim,input$kmax)
    p <- plot(roc1,print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
         grid.col=c("green", "red"), max.auc.polygon=TRUE,
         auc.polygon.col="skyblue", print.thres=TRUE)

    p
  })
  
  
  cluster_fun <- reactive({
    if(input$method == "LinkModule"){
      cluster1<-LinkModule(TF_miRNA(),TF_mRNA(),miRNA_mRNA(),mRNA_mRNA(),miRNA(),TF(),mRNA(),network_similarity())
    } else if(input$method == "LinkPartition"){
      cluster1<-LinkPartition(TF_miRNA(),TF_mRNA(),miRNA_mRNA(),mRNA_mRNA(),miRNA(),TF(),mRNA(),network_similarity())
    } else if(input$method == "DBSCAN"){
      cluster1<-dbscancluster(TF_miRNA(),TF_mRNA(),miRNA_mRNA(),mRNA_mRNA(),miRNA(),TF(),mRNA(),network_similarity())
    }
  })
  
  result_fun <- reactive({
    if(input$method == "LinkModule"){
      result <- read.table("cluster.result_LinkModule.txt",header = F, sep = "\t")
    } else if(input$method == "LinkPartition"){
      result <- read.table("cluster.result_LinkPartition.txt",header = F, sep = "\t")
    } else if(input$method == "DBSCAN"){
      result <- read.table("cluster.result_DBSCAN.txt",header = F, sep = "\t")
    }
    result
  })
  
  target_fun <- reactive({
    if(input$method == "LinkModule"){
      target <- read.table("targetgene.txt",header = F, sep = "\t")
    } else if(input$method == "LinkPartition"){
      target <- read.table("targetgene_LinkPartition.txt",header = F, sep = "\t")
    } else if(input$method == "DBSCAN"){
      target <- read.table("targetgene_DBSCAN.txt",header = F, sep = "\t")
    }
    target
  })
  
  go_fun <- reactive({
    go_terms <- goanalysis(target_gene(),input$method)
    go_terms
  })
  
  keggdown_fun <- reactive({
    keggdown <- read.csv("KEGG.csv",header = T,sep = ',')
    keggdown
  })
  
  #########output########
  
  output$net <- renderTable({
    dat = TF_miRNA()
  })
  
  output$ROC <- renderPlot({
    plot_fun_train()
  })
  
  observeEvent(input$cluster_button, {
    if(input$method == "LinkModule"){
      cluster1<-LinkModule(TF_miRNA(),TF_mRNA(),miRNA_mRNA(),mRNA_mRNA(),miRNA(),TF(),mRNA(),network_similarity())
    } else if(input$method == "LinkPartition"){
      cluster1<-LinkPartition(TF_miRNA(),TF_mRNA(),miRNA_mRNA(),mRNA_mRNA(),miRNA(),TF(),mRNA(),network_similarity())
    } else if(input$method == "DBSCAN"){
      cluster1<-dbscancluster(TF_miRNA(),TF_mRNA(),miRNA_mRNA(),mRNA_mRNA(),miRNA(),TF(),mRNA(),network_similarity())
    }
  })
  
  # go <- observeEvent(input$goterm_button, {
  #   go <-goanalysis()
  #   go
  # })
  # 
  output$summary <- renderDataTable({
    dat= go_fun()
    datatable(dat)
    # summary(dat)
  })
  
  output$download_result <- downloadHandler(
    filename = function(){
      paste0("clustering_result", ".txt")
    },
    contentType = "txt",
    content = function(file){
      dat1 = result_fun()
      fwrite(dat1,file,col.names = FALSE)
    }
  )
  
  output$download_target <- downloadHandler(
    filename = function(){
      paste0("targetgene", ".txt")
    },
    contentType = "txt",
    content = function(file){
      dat2 = target_fun()
      fwrite(dat2,file,col.names = FALSE)
    }
  )
  
  output$download_kegg <- downloadHandler(
    filename = function(){
      paste0("keggpathway", ".csv")
    },
    contentType = "csv",
    content = function(file){
      dat3 = keggdown_fun()
      write.csv(dat3, file)
    }
  )
  session$onSessionEnded(function() {
    stopApp()
  })
}
shinyApp(ui = ui, server = server)


library(shiny)
library(ggplot2)
library(colourpicker) # you might need to install this.
library(rsconnect)
library(xtable)
library(tidyverse)
library('RColorBrewer')
library(DT)
library(shinyWidgets)
library(ggbeeswarm)

options(shiny.maxRequestSize=40*1024^2) 
ui <- fluidPage(
  mainPanel(
    tabsetPanel(
      tabPanel("Sample Information Exploration",
        sidebarPanel(fileInput("input_1", paste0("Load sample data"),accept = c(".tsv",".csv")),
                     radioButtons("histo", "Chose a category to make a histogram of:", 
                                  c("PMI","age_death", "RIN","mRNA_seq_reads")),
                     colourInput("histcolor", "Histogram color"),
                     actionButton("go", "Go")
        ), #sidebarPanel
        mainPanel(
          tabsetPanel(
            tabPanel("Summary",
                     p("This is a summary of the uploaded data."),
                     tableOutput("sample_1")
            ),#tabPanel
            tabPanel("Sortable Data Table",
                     p("Click the arrows next to a column name to sort"),
                     DT::dataTableOutput("sample_2")
            ),#tabPanel
            tabPanel("Histogram",
                     p("Here's a histogram of the data of your chosing. Please give it a moment to load."),
                     plotOutput("sample_3")
            )#tabPanel
          )#tabsetPanel
        ), #mainPanel
      ), #Sample information tabPanel
      tabPanel("Counts Matrix Exploration",
               sidebarPanel(fileInput("input_2", paste0("Load sample data"),accept = c(".tsv",".csv")),
                            sliderInput("percentile_v", min = 0.11, max = 3.2, 
                                        label = "Include genes with at least X percentile of variance. Chose X: ",value = 1.2, step = 0.2),
                            sliderInput("percentile_z", min = 5, max = 60, 
                                        label = "Include genes with at least X samples that are non-zero. Chose X: ",value = 30, step = 5),
                            selectInput("x", "Which PC do you want on the x-axis", 
                                        choices = c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10'), selected = 1),
                            selectInput("y", "Which PC do you want on the y-axis", 
                                        choices = c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10'), selected = 1),
                            actionButton("go", "Go")
                            ),#sidebarPanel
               mainPanel(
                 tabsetPanel(
                   tabPanel("Summary of Filtering",
                            p("This is a summary of the uploaded data."),
                            tableOutput("counts_1")
                   ),#tabPanel
                   tabPanel("Scatter Plots",
                            p('This is a graph of median count vs variance. Please give it a moment to load (about 10 seconds).'),
                            plotOutput("counts_2"),
                            p('This is a graph of median count vs number of zeros. Please give it a moment to load (about 10 seconds).'),
                            plotOutput("counts_3")
                   ),#tabPanel
                   tabPanel("Heat Map",
                            p('Tab with a clustered heatmap of counts remaining after filtering. Please give it a moment to load (about 10 seconds).'),
                            plotOutput("counts_4")
                   ),#tabPanel
                   tabPanel("PCA Plot",
                            p('This is a scatter plot of principal component analysis projections. Please give it a moment to load (about 10 seconds).'),
                            plotOutput("counts_5")
                   )#tabPanel
                 )#tabsetPanel
               ),#mainPanel
      ),#Counts Matrix tabPanel
      tabPanel("Differential Expression",
               sidebarPanel(fileInput("input_3", paste0("Load sample data"),accept = c(".tsv",".csv")),
                            p('A volcano plot can be generated with "log2 fold-change" on the x-axis and "p-adjusted" on the y-axis.'),
                            radioButtons("x_name", "Choose the column for the x-axis", 
                                         c("baseMean","HD.mean","Control.mean","log2FoldChange","lfcSE","stat","pvalue","padj"),
                                         selected = "log2FoldChange"),
                            radioButtons("y_name", "Choose the column for the y-axis", 
                                         c("baseMean","HD.mean","Control.mean","log2FoldChange","lfcSE","stat","pvalue","padj"),
                                         selected = "padj"),
                            colourInput("color1", "Base point color"),
                            colourInput("color2", "Highlight point color"),
                            sliderInput("magnitude", min = -25, max = 0, 
                                        label = "Select the magnitude of the p adjusted coloring:",value = -15, step = 5),
                            actionButton("go", "Go")
               ),#sidebar
               mainPanel(
                 tabsetPanel(
                   tabPanel("Sortable Table",
                            p('Click the arrows next to a column name to sort'),
                            DT::dataTableOutput("difexp_1")
                   ),#tabPanel
                   tabPanel("Volcano Plot",
                            p('This is a volcano plot.'),
                            plotOutput('difexp_2')
                   ),#tabPanel
                   tabPanel("Filtered Table",
                            tableOutput('difexp_3')
                   ), #tabPanel
                 ),#tabsetPanel
               ), #mainPanel
      ),#DESeq TabPanel
      tabPanel("Visualization of Individual Gene Expression",
               sidebarPanel(fileInput("input_4_sample", paste0("Load sample data csv"),accept = c(".tsv",".csv")),
                            fileInput("input_4_counts", paste0("Load normalized counts csv"),accept = c(".tsv",".csv")),
                            selectInput("categorical_value", label = ("Pick a category from the sample data."), 
                                        choices = c('Control vs. HD Samples' = 1, 'Age at Death (36-55, 55-66, 66-75, 75-106)' = 2)),
                            helpText("Input one of the genes found in the counts matrix"),
                            searchInput("inputid",
                                        label = NULL,
                                        value = "ENSG00000000003.10",
                                        placeholder = "ENSG00000000003.10",
                                        btnSearch = NULL,
                                        btnReset = NULL,
                                        resetValue = "",
                                        width = NULL
                            ),
                            radioButtons("plot_type","What kind of plot would you like?",c('Barplot'=1,'Boxplot'=2,'Violin plot'=3,'Beeswarm'=4)),
                            actionButton("go", "Go")
               ),#sidebar
               mainPanel(
                 tabsetPanel(
                   tabPanel("Samples",
                            DT::dataTableOutput('cyo_2')
                   ),#tabPanel
                   tabPanel("Normalized Counts Table",
                            p('Click the arrows next to a column name to sort'),
                            DT::dataTableOutput("cyo_1")
                   ),#tabPanel
                   tabPanel("Plot",
                            plotOutput("cyo_4")
                   ),#tabPanel
                 ),#tabsetPanel
               ), #mainPanel
      ),#Chose your own Adventure TabPanel
    ) #tabsetPanel
  )#mainpanel
  )#fluidPage




server <- function(input, output, session){
  #Sample Tab Functions: 
  load_data <- reactive({
    read.csv(input$input_1$datapath,row.names = 1) #Reactive to file upload
  })
  
  sample_summary<- function(dataf){
    categories <- colnames(dataf) #find all the column names 
    type <- c(typeof(dataf$sample_id),typeof(dataf$PMI),typeof(dataf$age_death),
              typeof(dataf$RIN),typeof(dataf$mRNA_seq_reads)) #find the type of data stored in each column
    mean <-c("identifier",round(mean(dataf$PMI %>% na.omit),3),round(mean(dataf$age_death),3),
             round(mean(dataf$RIN),3),round(mean(dataf$mRNA_seq_reads),3)) #find the mean of each numerical column rounding to 3 digits
    sd <- ifelse(sapply(dataf, is.numeric) == TRUE,round(apply(dataf,2,sd,na.rm=TRUE),3),"identifier") #repeat with standard deviation
    summary <- data.frame(categories,type,mean,sd) #compile it into one dataframe
    return(summary)
  }
  
  histogram_plot <- function(df,category,bins){
    df$category<- as.numeric(df[[category]]) #use the [[]] to change a string identifier into a column name 
    p <- ggplot(df, aes(category)) +
      geom_histogram(fill = "blue",bins = bins) #make a histogram of this data 
    return(p)
  }
  
  output$text1 <- renderText({paste("Please give it a little time to process.")}) #Just a helper text
  output$sample_1 <- renderTable({req(input$input_1)  #to create the summary table output
    df<- load_data() 
    p<- sample_summary(df)
    return(p)
  }, height = 300)
  
  output$sample_2 <- DT::renderDataTable({req(input$input_1) #to create a sortable dataframe
    DT::datatable(load_data())})
  
  output$sample_3 <- renderPlot({req(input$input_1) #to make the histogram
    df <- load_data()
    p<- hist(as.numeric(df[[input$histo]]), xlab= paste(input$histo), main = paste("Histogram of ",input$histo), col = input$histcolor,
             breaks=10) 
    return(p)},
    height = 300)
  
  #Counts Tab 

  
  load_data_counts <- reactive({
    read.csv(input$input_2$datapath,row.names = 1)
  })
  
  filter_counts <- function(counts,percentile_v,percentile_z){
    df_variance <- mutate(counts,variance = (apply(counts,1,sd)/apply(counts,1,mean))>percentile_v) #min = 0.0897, max = 3.12 filter(variance ==TRUE)%>%  select(-variance)
    df_nonzero <- mutate(df_variance,zeros = rowSums(df_variance == 0)<percentile_z) #min = 0.0897, max = 3.12 filter(zeros ==TRUE)%>% select(-zeros)
    return(df_nonzero)
  }
  
  
  summary_counts <- function (counts,filtered){
    filtered <- filtered %>% filter(variance ==TRUE) %>% select(-variance) %>% filter(zeros ==TRUE)%>% select(-zeros)
    options(scipen = 50)
    samples <- ncol(counts)
    genes <- nrow(counts)
    excluded<-genes - nrow(filtered)
    new_num_genes<- nrow(filtered)
    df<- data.frame("Original number of samples"=c(samples,''),"Original number of genes"=c(genes,''),"Genes not passing the filter"=c(excluded,excluded/genes),"Genes passing the filter"=c(new_num_genes,new_num_genes/genes))
    rownames(df)<-c('Numbers',"Percentage of Total Genes")
    return(df %>%  mutate_if(is.numeric, round, digits = 3))
  }
  
  
  
  #Counts Tab 2
  scatter_cv<- function(filtered){
    df<- mutate(filtered,median = apply(filtered,1,median)) %>%  
      mutate(variance_cv = (apply(filtered,1,sd)/apply(filtered,1,mean)))%>% mutate(Filtered = ifelse(variance == FALSE | zeros ==FALSE ,'Filtered Out','Not Filtered Out')) #select(median,variance) 
    p<-ggplot(df,aes(median,variance_cv))+geom_point(aes(color = factor(Filtered)),alpha = 0.5)+ scale_x_continuous(trans='log2')+
      scale_color_manual(values = c('Filtered Out' = "mediumpurple1", "Not Filtered Out" = "mediumpurple4"))+
      labs(title="Plot of median counts vs Variance", x ="Log scale median counts", y = "Variance")
    return(p)
  } 
  
  scatter_nz <-function(filtered,percentile_z){
    df_nonzero <- mutate(filtered, median = apply(filtered,1,median)) %>% mutate(zeros_nz = rowSums(filtered == 0)) %>% 
      mutate(Filtered = ifelse(variance == FALSE | zeros ==FALSE ,'Filtered Out','Not Filtered Out'))
    p<-ggplot(df_nonzero,aes(median,zeros_nz))+geom_point(aes(color = factor(Filtered)),alpha = 0.5)+scale_x_continuous(trans='log2')+
      scale_color_manual(values = c('Filtered Out' = "palevioletred3", "Not Filtered Out" = "palevioletred4"))+
      labs(title="Plot of median counts vs Zeros", x ="Log scale median counts", y = "Number of Zeros ")
    return(p)
  }
  
  heatmap_plot<- function(f){
    f<- f%>%filter(variance ==TRUE) %>% select(-variance) %>% filter(zeros ==TRUE)%>% select(-zeros)  #filter only true values from filter function
    f[f==0] <- NA #replace all 0s with NA
    f<- f  %>%log10() %>% as.matrix() %>%  na.omit() %>% heatmap() #take the log, convert it to a matrix, omit NAs and make it a heatmap 
    return(p)
  }
  
  pca_plot<- function(counts,x,y){
    pca_results <- prcomp(scale(t(counts)), center=FALSE, scale=FALSE) #scale and transpose the counts
    matrix <- do.call(rbind.data.frame, pca_results) #binds the pca results 
    matrix_new <- mutate (matrix, col1 = rownames(matrix)) #turns rownames to column
    standard_deviation <- (pca_results$sdev) #saves the SD
    variance <- standard_deviation^2 #calculates variance
    total <- sum(variance) #finds total variance 
    variance_explained <- variance/total #divides each vairance by the table 
    new_tibble <- tibble(variance_explained = variance_explained) #making it a tibble 
    pc <- mutate(new_tibble,principal_components = seq.int(nrow(new_tibble)))
    pc$principal_components<-sub("^","PC",pc$principal_components) #names the PCs
    pc_cumulative <- mutate(pc,cumulative = cumsum(pc$variance_explained)) #cumulative sum of variance 
    x_plot <- matrix_new[,x] #x axis data
    y_plot <- matrix_new[,y] #y axis data
    x_variance <- round(pc_cumulative$variance_explained[strtoi(str_sub(x, start= -1))] *100,3) #turns x variance into a percentage
    y_variance <- round(pc_cumulative$variance_explained[strtoi(str_sub(y, start= -1))] *100,3) #turns y variance into a percentage 
    p <- ggplot()+geom_point(aes(x_plot,y_plot)) + labs(title = paste(x,"versus ", y), x= paste(x, ' (Variance Explained: ',x_variance,"%)"),
                                                        y = paste0(y,' (Variance Explainied: ',y_variance,"%)"))
    return(p)
  }
  
  output$counts_1 <- renderTable({req(input$input_2)
    df <- load_data_counts()
    filter <- filter_counts(df,input$percentile_v,input$percentile_z)
    p <-summary_counts(df,filter)
    return(p)},rownames = TRUE, height = 400)
  output$counts_2 <- renderPlot({req(input$input_2)
    df <- load_data_counts()
    filter <- filter_counts(df,input$percentile_v,input$percentile_z)
    p <- scatter_cv(filter) 
    return(p)},height = 400)
  output$counts_3 <- renderPlot({req(input$input_2)
    df <- load_data_counts()
    filter <- filter_counts(df,input$percentile_v,input$percentile_z)
    p <- scatter_nz(filter, input$percentile_z)
    return(p)}, height =400)
  output$counts_4 <- renderPlot({req(input$input_2)
    df <- load_data_counts()
    filter <- filter_counts(df,input$percentile_v,input$percentile_z)
    p<- heatmap_plot(filter)
    return(p)}, height = 400)
  output$counts_5 <- renderPlot({req(input$input_2)
    df <- load_data_counts()
    p <- pca_plot(df,input$x,input$y)
    return(p)}, height = 400)
  
  
  #DESeq Tab 
  
  load_data_difexp <- reactive({
    read.csv(input$input_3$datapath,row.names = 1)
  })
  
  draw_table <- function(dataf, slider) {
    dataf<-subset(dataf,padj< 1*10^slider) #subsets the dataf only padj passing the filter
    #is.num <- sapply(dataf, is.numeric)
    #dataf[is.num] <- lapply(dataf[is.num], round, 7) 
    dataf<- mutate(dataf,pvalue = format(pvalue, scientific = TRUE)) %>% mutate(padj = format(padj, scientific = TRUE)) #makes pvalue and padj into scientific notation
    return(dataf)
  }
  
  volcano_plot <- function(dataf, x_name, y_name, slider, color1, color2) {
      p <- ggplot(dataf, aes(x = !!sym(x_name),
                             y = -log10(!!sym(y_name)))) +
        geom_point(aes(color = !!sym(y_name) < 1 * 10 ^ (as.numeric(slider)))) +
        theme_bw() +
        scale_color_manual(values = c(color1, color2)) +
        theme(legend.position = "bottom") +
        labs(color = paste0(y_name, " < 1 × 10^", slider))
      return(p)
  }
  
  output$difexp_1 <- DT::renderDataTable({req(input$input_3)
    DT::datatable(load_data_difexp())})
  output$difexp_2 <- renderPlot({
    req(input$input_3)
    df <- load_data_difexp()
    p <-volcano_plot(load_data_difexp(),
                     input$x_name,
                     input$y_name,
                     input$magnitude,
                     input$color1,
                     input$color2)
    return(p)
  }, height = 700)
  
  output$difexp_3 <- renderTable({req(input$input_3)
    df <- load_data_difexp()
    p <-draw_table(df,input$magnitude)
    return(p)})
   
  #Chose your own Tab 
  load_counts_cyo <- reactive({
    read.csv(input$input_4_counts$datapath,row.names = 1)
  })
  
  load_sample_cyo <- reactive({
    read.csv(input$input_4_sample$datapath,row.names = 1)
  })
  
  gene<- reactive({
    input$inputid
  })
  
  
  age_division<- function(sample,norm,gene_id){
    dataf<-subset(norm,gene==gene_id) #the next 8 lines subsets the data into multiple dataframes depending on the age category
    age_1 <- filter(sample,age_death>=36) %>% filter(age_death<55)
    a1_columns <- age_1$sample_id
    age_2 <-filter(sample,age_death>=55) %>% filter(age_death<66)
    a2_columns <- age_2$sample_id
    age_3 <-filter(sample,age_death>=66) %>% filter(age_death<75)
    a3_columns <- age_3$sample_id
    age_4 <-filter(sample,age_death>=75) 
    a4_columns <- age_4$sample_id
    a1 <- dataf%>% select(all_of(a1_columns)) %>% unlist() %>% list() %>%lapply(function(x) {x[x!=0]}) #saves all the counts for each in a list, removing 0s 
    a2 <- dataf %>% select(all_of(a2_columns)) %>% unlist() %>% list() %>%lapply(function(x) {x[x!=0]})
    a3 <- dataf%>% select(all_of(a3_columns)) %>% unlist() %>% list() %>%lapply(function(x) {x[x!=0]})
    a4 <- dataf %>% select(all_of(a4_columns)) %>% unlist() %>% list() %>%lapply(function(x) {x[x!=0]})
    temp_1 <- data.frame(a1,letter = "Ages 36-55") #the first column is the count and the second column is the group
    colnames(temp_1) = c("counts","group")
    temp_2 <- data.frame(a2,letter = "Ages 55-66")
    colnames(temp_2) = c("counts","group")
    temp_3 <- data.frame(a3,letter = "Ages 66-75")
    colnames(temp_3) = c("counts","group")
    temp_4 <- data.frame(a4,letter = "Ages 75+")
    colnames(temp_4) = c("counts","group")
    ages <- rbind(temp_1,temp_2,temp_3,temp_4) #row bind the dataframes
    return(ages)
  }
  
  case_control<- function(sample,norm,gene_id){
    dataf<-subset(norm,gene==gene_id)
    c_vs_h <-mutate(sample,letter=substr(sample_id,1,1)) #whole function notes similar to one above
    c<- filter(c_vs_h,letter=='C') 
    c_columns <- c$sample_id
    h<- filter(c_vs_h,letter =="H")
    h_columns <- h$sample_id
    c_counts <- dataf %>% select(all_of(c_columns)) %>% unlist() %>% list() %>%lapply(function(x) {x[x!=0]})
    h_counts <- dataf %>% select(all_of(h_columns)) %>% unlist() %>% list() %>%lapply(function(x) {x[x!=0]})
    temp_h<-data.frame(h_counts,letter = "H")
    temp_c<- data.frame(c_counts,letter = "C")
    colnames(temp_h) = c("counts","group")
    colnames(temp_c) = c("counts","group")
    case_control_df <- rbind(temp_h,temp_c)
    return(case_control_df)
  }
  
  bar_plot <-function(dataframe){
    p<- ggplot(dataframe, aes(counts, fill = group)) + 
      geom_histogram(alpha = 0.5, bins = 30, aes(y = ..density..), position = 'identity')
    return(p)
    }
  
  bee_plot<-function(dataframe){
    p<-ggplot(dataframe, aes(x=group, y=counts, color = group)) + geom_beeswarm()+xlab("Cases vs. Controls")
    return(p)
  }
  
  box_plot <- function(dataframe){
    p <- ggplot(dataframe, aes(x=group, y=counts, fill=group)) + 
      geom_boxplot(alpha=0.3) +
      theme(legend.position="none")
    return(p)
  }
  violin_plot <- function(dataframe){
    p<- ggplot(dataframe, aes(x=group, y=counts, fill=group)) + geom_violin() +xlab("Cases vs. Controls")
    return(p)
  }
  
  

  
  output$cyo_1 <- DT::renderDataTable({req(input$input_4_counts)
    DT::datatable(load_counts_cyo())})
  
  output$cyo_2 <- DT::renderDataTable({req(input$input_4_sample)
    DT::datatable(load_sample_cyo())})
  
  output$cyo_3<- renderTable({req(input$input_4_counts)
                            df_counts<- load_counts_cyo() %>% rownames_to_column("gene")
                            df_samples<- load_sample_cyo()
                            gene_name <- input$inputid
                            p<- age_division(df_samples,df_counts,gene_name)
                            return(p)})
  output$cyo_4 <- renderPlot({req(input$input_4_counts)
                          df_counts<- load_counts_cyo() %>% rownames_to_column("gene")
                          df_samples<- load_sample_cyo()
                          gene_name <- gene()
                          if(input$categorical_value ==1){
                            data<- case_control(df_samples,df_counts,gene_name)
                          }
                          else{
                            data<- age_division(df_samples,df_counts,gene_name)
                          }
                          if(input$plot_type==1){
                            p<-bar_plot(data)
                          }
                          else if (input$plot_type==2){
                            p<-box_plot(data)
                          }
                          else if (input$plot_type==3){
                            p<- violin_plot(data)
                          }
                          else{
                            p<- bee_plot(data)
                          }
                          return(p)})

 
}


shinyApp(ui = ui, server = server)


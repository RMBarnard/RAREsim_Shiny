library(shiny)
library(RAREsim)
library(shinythemes)
library(DT)
library(dplyr)
library(ggplot2)

data("nvariant_afr")
data("afs_afr")

# Define UI for RAREsim
ui <- fluidPage(
    theme = shinytheme("cerulean"),
    navbarPage("RAREsim",
               id = "intabset",
        #Everything in this tab panel is for the execution of the main final product.
        tabPanel("Main",
            fluidRow(
                column(
                    width = 3,
                    wellPanel(
                        radioButtons("variant_parameters",
                                     "Nvariant Parameters",
                                   c("Default"= "default",
                                     "Custom" = "custom",
                                     "Fit" = "fit")),
                                    
                        conditionalPanel(condition = "input.variant_parameters == 'fit' && !output.nvar_uploaded",
                                         span(print("Warning: No input file detected for Nvariants.
                                                     Navigate to the 'Fit Target Data' tab to upload your data."), style="color:red")),
                                
                        radioButtons("afs_parameters_final",
                                     "AFS Parameters",
                                   c("Default"= "default",
                                     "Custom" = "custom",
                                     "Fit" = "fit")),
                                
                        conditionalPanel(condition = "input.afs_parameters_final == 'fit' && !output.afs_uploaded",
                                         span(print("Warning: No input file detected for AFS.
                                                     Navigate to the 'Fit Target Data' tab to upload your data."), style="color:red")),
                                
                        conditionalPanel(condition = "input.afs_parameters_final != 'fit'",
                                         radioButtons("bin_select_final",
                                                      "MAC Bins",
                                                    c("Default Bins" = "default",
                                                      "Custom Bins" = "custom"))),
                                
                        conditionalPanel(condition = "input.bin_select_final == 'custom' && input.afs_parameters_final != 'fit'",
                                         helpText("To enter custom bins, navigate to the AFS tab under 'Functions'")),
                        
                        conditionalPanel(condition = "input.variant_parameters == 'default' || input.afs_parameters_final == 'default' ",
                                         radioButtons("region",
                                                      "Region of Ancestry",
                                                     c("African" = "AFR",
                                                       "East Asian" = "EAS",
                                                       "Non-Finnish European" = "NFE",
                                                       "South Asian" = "SAS")))
                    ),
                    actionButton("gotohelp",
                                 "Help")
                ),
                column(
                    width = 3,
                    wellPanel(
                        #Numeric Input boxes for region size in kilobases and for sample size
                        numericInput("kbsize_total",
                                     "Region Size(Kb)",
                                     value = 19),
                                    
                        numericInput("sample_size_total",
                                     "Simulation Size (# of Individuals)",
                                     value = 8000),
                                
                        #Warning messages for if the desired sample size is too large or too small
                        conditionalPanel(condition = "input.sample_size_total > 125000",
                                         span(print("Warning: We currently do not recommend simulating
                                                    sample sizes over 125,000"), style="color:red")),
                                
                        conditionalPanel(condition= "input.sample_size_total<2000",
                                         span(print("Warning: To simulate <2000 individuals, use RAREsim to simulate 2000 individuals
                                                    and randomly downsample to the desired size"), style="color:red")),
                                
                        conditionalPanel(condition = "input.variant_parameters =='custom'",
                                         numericInput("phi_final",
                                                      HTML("Phi (&Phi;):"),
                                                      value = .6,
                                                      min = 0),
                                         conditionalPanel(condition = "input.phi_final <= 0",
                                                          span(print("Warning: Parameter out of range. Phi must be greater than 0."), style = "color:red")),
                                         numericInput("omega_final",
                                                      HTML("Omega (&omega;)"),
                                                      value = 0.45),
                                         conditionalPanel(condition = "input.omega_final >= 1 || input.omega_final <= 0",
                                                          span(print("Warning: Parameter out of range. Omega must be between 0 and 1"),style = "color:red"))
                        ),
                        conditionalPanel(condition = "input.afs_parameters_final =='custom'",
                                         numericInput("alpha_final",
                                                      HTML("Alpha (&alpha;)"),
                                                      value = NA),
                                         conditionalPanel(condition = "input.alpha_final <= 0",
                                                          span(print("Warning: Parameter out of range. Alpha must be greater than 0."), style = "color:red")),
                                         numericInput("beta_final",
                                                      HTML("beta (&beta;)"),
                                                      value = NA),
                                         numericInput("b_final",
                                                      HTML("<i>b</i>"),
                                                      value = NA))
                    )
                ),
                column(
                    width = 6,
                    wellPanel(
                        #The main panel was here
                        h3("Final Results"),
                        tableOutput("final_table"),
                        downloadButton("TXTDown", label = "Download Table")
                    ),
                    wellPanel(
                        textOutput("main_help")
                    )
                )
            )
        ),
        navbarMenu("Functions",
            #All work for NVariants is in this panel
            tabPanel("NVariants",
                # Application title
                titlePanel("NVariants"),
                # Sidebar with various inputs for nvariants function
                fluidRow(
                    column(
                        width = 3,
                        wellPanel(
                            #Radio Buttons for selection of either default or custom parameters
                            radioButtons("nvar_parameters",
                                         "Parameters",
                                       c("Default" = "default",
                                         "Custom" = "custom")),
                                               
                            #Radio Button panel that will appear when default parameters are selected
                            #to select an ancestry population
                            conditionalPanel(condition = "input.nvar_parameters == 'default'",
                                             radioButtons("nvar_region",
                                                          "Region of Ancestry",
                                                        c("African" = "AFR",
                                                          "East Asian" = "EAS",
                                                          "Non-Finnish European" = "NFE",
                                                          "South Asian" = "SAS"))),
                                               
                            #Numeric Input boxes for phi and omega 
                            #that will appear if custom parameters are selected
                            conditionalPanel(condition = "input.nvar_parameters =='custom'",
                                             numericInput("phi",
                                                          HTML("Phi (&Phi;):"),
                                                          value = 0.6),
                                             conditionalPanel(condition = "input.phi <= 0",
                                                              span(print("Warning: Parameter out of range. Phi must be greater than 0."))),
                                             numericInput("omega",
                                                          HTML("Omega (&omega;)"),
                                                          value = 0.45),
                                             conditionalPanel(condition = "input.omega >= 1 || input.omega <= 0",
                                                              span(print("Warning: Parameter out of range. Omega must be between 0 and 1"),style = "color:red"))),
                            
                            #Numeric Input boxes for region size in kilobases and for sample size
                            numericInput("kbsize",
                                         "Region Size(Kb)",
                                         value = 19),
                            numericInput("nvar_size",
                                         "Sample Size",
                                         value = 8000),
                                               
                            #Warning messages for if the desired sample size is too large or too small
                            conditionalPanel(condition = "input.nvar_size > 125000",
                                             span(print("Warning: We currently do not recommend simulating
                                                         sample sizes over 125,000"), style="color:red")),
                            conditionalPanel(condition = "input.nvar_size<2000",
                                             span(print("Warning: To simulate <2000 individuals, use RAREsim to simulate 2000 individuals
                                                         and randomly downsample to the desired size"), style="color:red"))
                        )
                    ),
                    column(
                        width = 5,
                        wellPanel(
                            #Show a plot of the generated distribution
                            h3("Results"),
                            textOutput("variants")
                        ),
                        wellPanel(
                            textOutput("nvar_help")
                        )
                    )
                )
            ),
            #All work for AFS is in this panel
            tabPanel("AFS",
                titlePanel("Allele Frequency Spectrum"),
                fluidRow(
                    column(
                        width = 3,
                        wellPanel(
                            #Radio Buttons to allow for the selection of custom or default bins
                            radioButtons("bin_select",
                                         "Select Bins",
                                       c("Default Bins" = "default",
                                         "Custom Bins" = "custom")),
                            #This conditional panel will output an interactive table for the entry
                            #of custom bins if desired
                            conditionalPanel(condition = "input.bin_select == 'custom'",
                                             numericInput("bins",
                                                          "Number of MAC Bins",
                                                          value = 1),
                                                          DTOutput("my_datatable")),
                            #Radio Buttons to allow for the selection of default or custom parameters
                            radioButtons("afs_parameters",
                                         "Parameters",
                                       c("Default" = "default",
                                         "Custom" = "custom")),
                            #This conditional panel holds radio buttons to select an ancestry if
                            #default parameters are desired
                            conditionalPanel(condition = "input.afs_parameters == 'default'",
                                             radioButtons("afs_region",
                                                          "Region of Ancestry",
                                                        c("African" = "AFR",
                                                          "East Asian" = "EAS",
                                                          "Non-Finnish European" = "NFE",
                                                          "South Asian" = "SAS"))),
                            #This conditional panel holds a numeric input box to get the sample
                            #size for when default bins are selected. Warning messages are included
                            conditionalPanel(condition = "input.bin_select == 'default'", 
                                             numericInput("afs_size",
                                                          "Sample Size",
                                                          value = 8000),
                            conditionalPanel(condition = "input.afs_size < 2000",
                                             span(print("Warning: We currently do not recommend simulating
                                                         sample sizes less than 2000"), style="color:red"))),
                                          
                            #This conditional panel holds various numeric input boxes to
                            #allow for the input of custom parameters
                            conditionalPanel(condition = "input.afs_parameters =='custom'",
                                             numericInput("alpha",
                                                          HTML("Alpha (&alpha;)"),
                                                          value = 0),
                                             conditionalPanel(condition = "input.alpha <= 0",
                                                              span(print("Warning: Parameter out of range. Alpha must be greater than 0."), style = "color:red")),
                                             numericInput("beta",
                                                          HTML("beta (&beta;)"),
                                                          value = 0),
                                             numericInput("b",
                                                          HTML("<i>b</i>"),
                                                          value = 0))
                        )    
                    ),
                    column(width = 5,
                        wellPanel(
                            #The main panel will simply output the AFS table with proportions
                            h3("Results"),
                            tableOutput("afs_table")
                        ),
                        wellPanel(
                            textOutput("afs_help")
                        )
                    )
                )
            )
        ),
        navbarMenu("Fit Target Data",
            #All work for the fit_nvariants function is in this tab
            tabPanel("Fit_Nvar",
                titlePanel("Fit_Nvariants"),
                fluidRow(
                    column(
                        width = 3,
                        wellPanel(
                            fileInput("nvariant_input",
                                      "Upload TXT File",
                                      accept = ".txt",
                                      multiple = FALSE,
                                      buttonLabel = "Choose File to Upload"),
                            helpText("To use the Fit_Nvariant function, you will need
                                      to upload a txt file in the same format and with
                                      the same column names as the example table below"),
                            conditionalPanel("!output.nvar_uploaded",
                                             h4("Example Input"),
                                             tableOutput("nvariant_example")),
                            conditionalPanel("output.nvar_uploaded",
                                             h4("Uploaded Input"),
                                             tableOutput("uploaded_nvar"))
                        )
                    ),
                    column(
                        width = 9,
                        wellPanel(
                            h3("Output"),
                            verbatimTextOutput("fit_nvariant_results"),
                            plotOutput("nvariant_plot")
                        ),
                        wellPanel(
                            textOutput("fit_nvar_help")
                        )
                    )
                )
            ),
            #All work for the fit_AFS function is in this tab
            tabPanel("Fit_AFS",
                titlePanel("Fit_AFS"),
                fluidRow(
                    column(
                        width = 3,
                        wellPanel(
                            #File input that will accept a TXT file. Help notes
                            #are also included to make sure that the necessary file
                            #format is clear to users
                            fileInput("afs_input",
                                      "Upload TXT File",
                                      accept = ".txt",
                                      multiple = FALSE,
                                      buttonLabel = "Select File..."),
                            helpText("To use the Fit_AFS function you will need
                                      to upload a txt file in the same format and with
                                      the same column names as the example table
                                      below"),
                            #This is a table output with an example of what the
                            #input data should look like
                            conditionalPanel(condition = "!output.afs_uploaded",
                                             h4("Example Input"),
                                             tableOutput("afs_example")),
                            conditionalPanel(condition = "output.afs_uploaded",
                                             radioButtons("maf_bins",
                                                          "Were any of these bin sizes calculated based on MAF?",
                                                        c("Yes" = "yes",
                                                          "No" = "no")),
                                conditionalPanel(condition = "input.maf_bins == 'yes'",
                                                 numericInput("breakpoint",
                                                              "What is the lower bound of the first bin to use MAF calculations",
                                                              value = 0),
                                                 numericInput("target_size",
                                                              "How large was the sample size that your target data was collected from?",
                                                              value = 10000),
                                                              h4("Uploaded Input"),
                                                              tableOutput("uploaded_afs"))
                            )
                        )
                    ),
                    column(
                        width = 9,
                        wellPanel(
                            #Outputting the fitted parameters with a plot
                            h3("Output"),
                            verbatimTextOutput("fit_afs_results"),
                            plotOutput("afs_plot")
                        ),
                        wellPanel(
                            textOutput("fit_afs_help")
                        )
                    )
                )
            )
        ),
        tabPanel("ReadMe",
                 value = "readme",
            fluidRow(
                column(7,
                    wellPanel(
                        h3("ReadME"),
                        br(),
                        print("The RAREsim Shiny App can be used to create the expected rare 
                              allele count bins to be used in RAREsim simulations. 
                              Function parameters can be estimated from target data
                              within the 'Fit Target Data' tab. For fully detailed instructions click"),
                        actionLink("help",
                                   "here.",
                                   onclick ="window.open('https://github.com/RMBarnard/RAREsim_Shiny/blob/main/Shiny%20Instructions.pdf', '_blank')")
                    ),
                    img(src = "Instructions4.png", width = "99%", height = "99%", align = "left"),
                    wellPanel(
                        br(),
                        print("*If custom or default parameters are selected and you wish to use custom MAC bins, the bins can be input
                              on the AFS page under the 'Functions' tab")
                    )
                ),
                column(5,
                    wellPanel(
                        h3("Acknowledgements"),
                        HTML("<b>RAREsim was developed as a collaborative effort by:</b></br>
                             Megan Null, Jos&eacute;e Dupuis, Pezhman Sheinidashtegol, Ryan M Layer, ChristopherR. Gignoux, and Audrey E. Hendricks
                             </br></br><b>The RAREsim Shiny App was developed as a collaborative effort by:</b></br>
                             Ryan Barnard, Nour Ibrahim, Megan Null, and Jessica Murphy</br></br>
                             RAREsim development was supported by the National Human Genome Research Institute (R35HG011293 and U01HG009080 to AEH and CGR; U01GH009080-05S1 to CGR).
                             The RAREsim Shiny App development was supported by an Institutional Development Award (IDeA) from the National Institute of General Medical Sciences
                             of the National Institutes of Health under Grant #P20GM103408.")
                    ),
                    img(src='coi_logo.png',width = "75%",height = "75%", align = "center"),
                    img(src = 'CU_Denver_logo.png',width = "85%",height = "85%",  align = "center")
                )
            )
        )
    )
)

#Define the server logic that allows us to make our functions react to
#the given inputs
server <- function(input, output, session) {
    
    #Below is the work for the output of the nvariants function
    #We start by adding a "variants" object to our output list
    #We then have "if" statements to allow us to run the nvariants function
    #with whatever parameters are given. We then store the return value of the
    #nvariants function multiplied by the size of the region and store it as "x"
    output$variants <- renderText({
        #If default parameters are selected, we run the nvariant function based on ancestry region
        if(input$nvar_parameters == "default"){
            x <- input$kbsize*nvariant(N = input$nvar_size, pop = input$nvar_region)
        }
        #Otherwise, if custom parameters are selected, we run the nvariant function using the
        #given phi and omega values.
        else if(input$nvar_parameters == "custom"){
            x <- input$kbsize*nvariant(N = input$nvar_size, phi = input$phi, omega = input$omega)
        }
        #Here we are formatting y to be the same as x, just rounded to 2 decimal places
        y <- formatC(x, digits = 2, format = "f")
        paste("Total number of expected variants:", y)
        
    }) 
    
    #This is the main output of the AFS function. afs_table is a table with bins and proportions
    #that correspond with each bin
    output$afs_table <- renderTable({
        #If default bins are selected we will use the recommended default bins depending on sample size
        if(input$bin_select == "default"){
            if(input$afs_size >= 3500){
                mac <- data.frame(Lower = c(1, 2, 3, 6, 11, 21, floor(input$afs_size*2*.005) + 1),
                                  Upper = c(1, 2, 5, 10, 20, floor(input$afs_size*2*.005), floor(input$afs_size*2*.01))) 
            }
            else if(input$afs_size < 3500){
                mac <- data.frame(Lower = c(1, 2, 3, 6, floor(input$afs_size*2*.0025) + 1, floor(input$afs_size*2*.005) + 1),
                                  Upper = c(1, 2, 5, floor(input$afs_size*2*.0025), floor(input$afs_size*2*.005), floor(input$afs_size*2*.01)))
            }
            #We also need to check for which parameters to use based on the selected parameter types
            if(input$afs_parameters == "default"){
                temp <- afs(pop = input$afs_region, mac_bins = mac) 
            }
            else if(input$afs_parameters == "custom"){
                temp <- afs(alpha = input$alpha, beta = input$beta, b = input$b, mac_bins = mac)
            }
        }
        #This else statement contains the logic for when custom bins are selected
        else{
            mac <- v$data[1:input$bins,]
            if(input$afs_parameters == "default"){
                temp <- afs(pop = input$afs_region, mac_bins = mac) 
            }
            else if(input$afs_parameters == "custom"){
                temp <- afs(alpha = input$alpha, beta = input$beta, b = input$b, mac_bins = mac)
            }
        }
        #Here we format all of our columns properly before outputting our table
        temp[, 1:2] = format(round(temp[,1:2],0))
        temp[, 'Prop'] = format(round(temp[,'Prop'],4),nsmall =4)
        expr = temp
    })
    
    #V is a backing datafram that stores the values for the interactive table that allows
    #for the selection of custom bins for the AFS function
    v <- reactiveValues(data = { 
        data.frame(Lower = numeric(0),Upper = numeric(0)) %>% 
            add_row(Lower = rep(1,100),Upper = rep(1,100))
    })
    
    #output the datatable based on the dataframe (and make it editable)
    output$my_datatable <- renderDT({
        DT::datatable(v$data[1:input$bins,],editable = TRUE, options = list(dom = 't',pageLength = 100 , ordering=F))
    })
    
    observeEvent(input$my_datatable_cell_edit, {
        #Get all info from the edited cell
        info = input$my_datatable_cell_edit
        i = as.numeric(info$row)
        j = as.numeric(info$col)
        k = as.numeric(info$value)
        
        #Write values to backing data frame
        v$data[i,j] <- k
    })
    
    #afs_example is a simple table that will display an example of how an uploaded
    #dataset should look
    output$afs_example <- renderTable({
        temp <- afs_afr
        temp[, 1:2] = format(round(temp[,1:2],0))
        temp[, 'Prop'] = format(round(temp[,'Prop'],4),nsmall =4)
        expr = temp
    })
    
    #fit_afs_results calculates and prints the fitted b, alpha, and beta values from
    #the given target data
    output$fit_afs_results <- renderText({
        if(is.null(input$afs_input)){
            print("No file found")
        }
        else{
            df <- read.table(input$afs_input$datapath, sep = '\t', header = TRUE)
            af <- fit_afs(df)
            b <- formatC(af$b, digits = 4, format = "f")
            beta <- formatC(af$beta, digits = 4, format = "f")
            alpha <- formatC(af$alpha, digits = 4, format = "f")
            paste("b:", b,
                  "\nbeta:", beta,
                  "\nAlpha:", alpha)
        }
    })
    
    output$uploaded_afs <- renderTable({
        df <- read.table(input$afs_input$datapath, sep = '\t', header = TRUE)
        expr = df
    })
    
    output$uploaded_nvar <- renderTable({
        df <- read.table(input$nvariant_input$datapath, sep = '\t', header = TRUE)
        expr = df
    })
    
    #afs_plot is a plot that will show all of the target values and all of our fitted
    #values on the same plot to allow for visualization of how well the fitted
    #function matches the given target data
    output$afs_plot <- renderPlot({
        #This if statement is included to prevent unwanted errors from being displayed
        #before a file is uploaded.
        if(!is.null(input$afs_input)){
            df <- read.table(input$afs_input$datapath, sep = '\t', header = TRUE)
            af <- fit_afs(Observed_bin_props = df)
            df$Data <- 'Target Data'
            af$Fitted_results$Data <- 'Fitted Function'
            af_all <- rbind(df, af$Fitted_results)
            af_all$Data <- as.factor(af_all$Data)
            p <- ggplot(data = af_all)+
                geom_point(mapping= aes(x= Lower, y= Prop, col= Data), size = 3,alpha = .75)+
                labs(title = "Target Data vs. Fitted Points",
                     x= 'Lower Limit of MAC Bin',
                     y= 'Proportion of Variants')+
                theme(axis.title = element_text(size = 16),
                      axis.text = element_text(size = 12),
                      plot.title = element_text(size = 20),
                      plot.background = element_rect(color = "#d9d9d9",fill = "#f5f5f5",size = 1, linetype = 1),
                      panel.background = element_rect(fill = "#ffffff"),
                      panel.grid = element_line(color = "#f5f5f5"),
                      legend.background = element_rect(fill = "#ffffff"),
                      legend.key = element_rect(fill = "#ffffff"),
                      legend.title = element_text(size = 16),
                      legend.text = element_text(size = 14)
                )
            plot(p)
        }
    })
    
    output$nvariant_example <- renderTable({
        expr = nvariant_afr
    }) 
    
    output$fit_nvariant_results <- renderText({
        if(is.null(input$nvariant_input)){
            print("No file found")
        }
        else{
            df <- read.table(input$nvariant_input$datapath, sep = '\t', header = TRUE)
            af <- fit_nvariant(Observed_variants_per_kb = df)
            phi <- formatC(af$phi, digits = 4, format = "f")
            omega <- formatC(af$omega, digits = 4, format = "f")
            paste("Phi:", phi,
                  "\nOmega:", omega)
        }
    })
    
    output$nvariant_plot <- renderPlot({
        if(!is.null(input$nvariant_input)){
            df <- read.table(input$nvariant_input$datapath, sep  = '\t', header = TRUE)
            af <- fit_nvariant(df)
            p <- ggplot(data = df)+
                geom_line(mapping = aes(x = n, y = (af$phi * (n^af$omega)), col= 'red'),size = 1)+
                geom_point(mapping = aes(x = n, y = per_kb),size = 2.5)+
                labs(x = "Sample Size", y = "Variants per Kb",
                     title = "Given Points with Fitted Function")+ 
                theme(legend.position = "none",
                      axis.title = element_text(size = 16),
                      axis.text = element_text(size = 12),
                      plot.title = element_text(size = 20),
                      plot.background = element_rect(color = "#d9d9d9",fill = "#f5f5f5",size = 1, linetype = 1),
                      panel.background = element_rect(fill = "#ffffff"),
                      panel.grid = element_blank()
                )
            
            plot(p)   
        }
        
    })
    
    getDataNvar <- reactive({
        if(is.null(input$nvariant_input)){
            return(NULL)
        }
        else{
            return(TRUE)
        }
    })
    
    getDataAFS <- reactive({
        if(is.null(input$afs_input)){
            return(NULL)
        }
        else{
            return(TRUE) 
        }
    })
    
    output$nvar_uploaded <- reactive({
        return(!is.null(getDataNvar()))
    })
    
    output$afs_uploaded <- reactive({
        return(!is.null(getDataAFS()))
    })
    
    outputOptions(output, 'nvar_uploaded', suspendWhenHidden = FALSE)
    outputOptions(output, 'afs_uploaded', suspendWhenHidden = FALSE)
    
    data <- reactive({
        #Calculating number of variants
        if(input$variant_parameters == "default"){
            z <- input$kbsize_total*nvariant(N = input$sample_size_total, pop = input$region)
        }
        else if(input$variant_parameters == "custom"){
            z <- input$kbsize_total*nvariant(N = input$sample_size_total, phi = input$phi_final, omega = input$omega_final)
        }
        else if(input$variant_parameters == "fit"){
            df <- read.table(input$nvariant_input$datapath, sep = '\t', header = TRUE)
            af <- fit_nvariant(Observed_variants_per_kb = df)
            z <- input$kbsize_total*nvariant(N = input$sample_size_total, phi = af$phi, omega = af$omega)
        }
        
        
        #Calculating Bin Proportions
        if(input$afs_parameters_final == "fit"){
            df <- read.table(input$afs_input$datapath, sep = '\t', header = TRUE)
            af <- fit_afs(df)
            if(input$maf_bins == "yes"){
                a = input$breakpoint
                temp <- df[which(df$Lower >= a),1:2]
                temp[,2] <- floor((temp[,2]/input$target_size) * input$sample_size_total)
                temp[2:nrow(temp),1] <- temp[1:nrow(temp)-1,2]+1
                mac <- merge(x = df[which(df$Lower < a),1:2], y = temp, by = c("Lower","Upper"), all = TRUE)
            }
            else if(input$maf_bins == "no")
            {
                mac = df[,1:2]
            }
            
            y <- afs(b = af$b, alpha = af$alpha, beta = af$beta, mac = mac)
        }
        
        if(input$bin_select_final == "default"){
            if(input$afs_size >= 3500){
                mac <- data.frame(Lower = c(1, 2, 3, 6, 11, 21, floor(input$sample_size_total*2*.005) + 1),
                                  Upper = c(1, 2, 5, 10, 20, floor(input$sample_size_total*2*.005), floor(input$sample_size_total*2*.01))) 
            }
            else if(input$afs_size < 3500){
                mac <- data.frame(Lower = c(1, 2, 3, 6, floor(input$sample_size_total*2*.0025) + 1, floor(input$sample_size_total*2*.005) + 1),
                                  Upper = c(1, 2, 5, floor(input$sample_size_total*2*.0025), floor(input$sample_size_total*2*.005), floor(input$sample_size_total*2*.01)))
            } 
        }
        
        else if(input$bin_select_final == "custom"){
            mac <- v$data[1:input$bins,]
        }
        
        if(input$afs_parameters_final == "default"){
            y <- afs(mac_bins = mac, pop = input$region)
        }
        
        else if(input$afs_parameters_final == "custom"){
            y <- afs(mac_bins = mac, alpha = input$alpha_final, beta = input$beta_final, b = input$b_final)
        }
        
        
        #Creating Final Table
        x = data.frame(expected_variants(Total_num_var = z, mac_bin_prop = y))
    })
    
    output$final_table <- renderTable({data()})
    
    output$TXTDown <- downloadHandler(
        filename = function() {
            paste('RAREsim-', Sys.Date(), '.txt', sep = '')
        },
        content = function(file) {
            write.table(data(), file, sep = "\t", row.names = FALSE)
        }
    )
    output$main_help <- renderText({
        print("In the above table you can see the expected number of variants for your simulation size in each bin.
               The dowload button will give you this table as a .txt file to upload for simulation.")
    })
    
    output$afs_help <- renderText({
        print("In this table you can see the proportion of the total variants that we expect to see in each bin.")
    })
    
    output$nvar_help <- renderText({
        print("The number above is the total number of rare variants that we would expect to see over the given kilobase range and sample size.")
    })
    
    output$fit_nvar_help <- renderText({
        print("The plot above shows the given datapoints as black dots along with the function that we fitted from the given datapoints.")
    })
    
    output$fit_afs_help <- renderText({
        print("The plot above shows the given datapoints as orange dots, along with the points that our function generated from the target data shown as blue dots.")
    })
    
    observeEvent(input$gotohelp, {
        updateTabsetPanel(session, inputId = "intabset", selected = "readme")
    })
    
}
# Run the application 
shinyApp(ui = ui, server = server)

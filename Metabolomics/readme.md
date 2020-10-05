# A collection of R snippets for metabolomics data processing

## Table of Contents

* [Data processing](#data-processing)
   * [Split long metabolite names](#Split-long-metabolite-names)
* [Handling msp-files](#handling-msp-files)
   * [Read msp-files](#read-msp-files)
   * [Plot Spectra from msp-data style](#Plot-Spectra-from-msp-data-style)
   * [msp-Viewer](#msp-viewer)

## Data processing

This section reports some potential tools to process the data in R.

### Split long metabolite names
Usually you don't want to work with the full database ID but with a shorter metabolite name. To do this you can shorten the obtained metabolite identification at the first or better the second underscore, to keep the derivatisation status

#### Generate Demodata
```R

  library(dplyr)
  library(tidyr)
  ## Let's say we have a dataframe with some measurements:
  data = tibble(Compound = c("Pyruvic acid_(1MEOX)(1TMS)_MP_RI:1053_IDENT:A+D",
                             "Lactic acid_(2TMS)_MP_RI:1065_IDENT:B+C",
                             "Valine_(1TMS)_BP_RI:1087_IDENT:B+D",
                             "Alanine_(2TMS)_BP_RI:1106_IDENT:B+C",
                             "Glycine_(2TMS)_BP_RI:1124_IDENT:A+C",
                             "Leucine(1TMS)_BP_RI:1158_IDENT:B+C"),   # here we introduced an error by accident
                 Intensity = c(250000, 12045969, 104949, 56554, 3996205, 482383))

```

#### Option 1
```R
  ## OPTION 1 - Just take the name as entry before the second underscore  - Thanks to Tobi
  data_sep1A = mutate(data, Name = sub("(.*?_.*?)_(.*)","\\1", Compound))

  # doing the same using str_extract() in the stringr library (tidyverse)
  library(stringr)
  data_sep1B = mutate(data, Name = str_extract(Compound, ".*?_.*?(?=_.*)"))
```R

#### Option 2
```R
  ## OPTION 2 - Just split into all entities along the underscore
  # this might lead to problems, when there are missing or too much underscores, see for example Leucine
   data_sep2 = separate(data, Compound, into = c("Name", "Derivative",
                                                "Product", "RI", "Ident"), sep = "_") 

  # 2B - splitted data can then merged together again
  data_united = unite(data_sep2, Metabolite, c("Name", "Derivative", "Product"), sep = "_") 
```
#### Option 3
```R
  ## OPTION 3. Split directly at the second underscore 
   TO DO
```
#### Option 4
```R
  ## OPTION 4 - FAIL_SAFE APPROACH: Split before the RI:
  # Clearly this only works when "_RI:" is there   
   data_sep4 = separate(data, Compound, into = c("Name", "Details"), sep = "_(?=RI:)") 
```
[back to top](#table-of-contents)

## Handling msp-files
This section reports some ideas to work with mass spectra (msp) files.  

### Read msp-files
This function allows to import msp-files into R for further processing. It reads msp files and generates a dataframe. As the structure can be very different the names are extracted and put in the first column, the spectra in the third and the comments (everything in between) into the second column.  

```R
#
#  Function to simply import mass spectra files and converts these in a dataframe. 
#  It is designed to read a variety of different styles, so it doesn't require (a lot of) fixed entries. 
#  The idea is to get an overview of the file, to compare the informations, which is difficult otherwise.
#  The dataframe can be further modified or exported as file (e.g. csv or excel)
#  required library: readr, stringr, magrittr (or simply tidyverse) - I assumed you may need theses libraries anyway.
#
#  usage: target_df = read_msp(Inputfile) , where "Inputfile" is the full name (and maybe the path) to your .msp file
#         and target_df is the name of the dataframe you want to store the data in


read_msp = function(inputfile) {

## Read in the msp
input = read_file(inputfile) %>% 
        str_replace_all("\r*\n", "\n") %>% 
        str_remove("\n\n$")                             # remove empty line at the end
  
## transform the file for better readability:
# split the different entries, typically separated by a blank line
input2 = str_split(input, "\n\n", simplify = FALSE)  # generates list
input3 = unlist(input2)                                  # back to character vector

## extract the important information
# extract compound names: everything between "Name:" and the next line break
names = str_extract(input3, regex("(?<=^(Name|NAME|name): ).*(?=\n)", # somehow ignore case doesnt work?
                    ignore_case = TRUE))  

# the next steps aren't very elegant but work
# remove the names from the input (as this nicely ends with a line break)
# line breaks mess up with the extraction of the other features, 
# so replace line-breaks by "___" to restore them later 

input_wo_names = str_remove(input3, regex("^(Name|NAME|name): .*\n", ignore_case = TRUE)) %>%  
                            str_replace_all("\n","___")

# extract comments: everything from after the name and (including "Num Peaks: NNN)
# then restore line breaks
comments = str_replace_all(str_extract(input_wo_names,
                                       regex(
                                       ".*Num Peaks: *\\d{1,3}", ignore_case = TRUE)), 
                               "___", "\n")

#extract peaks: everything after "Num Peaks: NNN" to the end
peaks = str_extract(input_wo_names,
                    regex("(?<=Num Peaks: \\d{1,3}___).*" , ignore_case = TRUE))
# clean up, remove unneeded symbols
peaks = str_replace_all(peaks,
                        c("___" = "",     # old line breaks
                          "; *" = ";",    # empty space after semicolon
                          "^ *" = "",     # empty space in front
                          ";$" = ""))     # ending semicolon

# create a dataframe by combining the 3 vectors:
msp_as_df = tibble("Name" = names,
                   "Comments" = comments,
                   "Peaks" = peaks) 

# check the results:
if (is.data.frame(msp_as_df) & all(colnames(msp_as_df) == c("Name", "Comments", "Peaks")) == TRUE) {
  print(paste0("Done! ", length(unique(msp_as_df$Name)),
               " unique compounds with ", length(msp_as_df$Peaks),
               " spectra imported"))
  
  return(msp_as_df)
  } else {
   print("Error, conversion failed!")     # doesn't really work, even when something wnet wrong a dataframe (with NAs) is generated
   return(msp_as_df) }
}
```
[back to top](#table-of-contents)


### Plot Spectra from msp-data style

```R
	# function: Plots a spectrum from msp format.
# author: Tobias Opialla
# required packages: ggplot2
# expects: string in msp-Format (m/z and intensity separated by space(s), pairs separated by semicolon - e.g.: "60   27;  61   19;  62    1;  65    1;  66    8")
# Options:
#    mytitle              - "name of the compound shown in the title " (character)
#    print_plot           - should plot be shown? (TRUE / FALSE)
#    return_plot          - should plot be exported as ggplot object? (TRUE / FALSE)
#    return_dataframe     - should the underlying dataframe be exported as struvtured table? (TRUE / FALSE)
#    normalize_Intensity  -
#    upper_limit_mzrange  - allows to set a limit for the x-axis (numeric)

plotspectra=function(stringfromMSP, mytitle="substance", print_plot=T, return_plot=F, 
                     return_dataframe=F, normalize_Intensity=TRUE, upper_limit_mzrange=0){
 #organize error handling
 if(return_plot & return_dataframe){
    print("please decide if you want either the data frame or the plot returned, the plot can be printed anyways")
    return(NULL)
  }

  if(grepl(';',stringfromMSP)){
    splitted=unlist(strsplit(stringfromMSP,';'))
    splitted=gsub("\\n",'',splitted)
    splitted=sub("^\\s+", "", splitted)

    df=data.frame(m.z=gsub("(\\d+)\\s+(\\d+)",'\\1',splitted),
                  intensity=gsub("(\\d+)\\s+(\\d+)",'\\2',splitted),stringsAsFactors = F)
    df=data.frame(lapply(df,as.numeric))
  }else{
    splitted=unlist(strsplit(stringfromMSP,' '))
    df=data.frame(m.z=gsub("(\\d+.0):(\\d+)",'\\1',splitted),
                  intensity=gsub("(\\d+.0):(\\d+)",'\\2',splitted),stringsAsFactors = F)
    df=data.frame(lapply(df,as.numeric))
    if(normalize_Intensity){
      df$intensity=(df$intensity/(max(df$intensity)))*1000
      df=df[df$intensity>5,]
    }
  }
  #filter empty lines
  df=df[!is.na(df$m.z),]

 #create plot
  p=ggplot(df,aes(x=m.z,y=intensity))+
    geom_bar(stat="identity")+
    ggtitle(mytitle)
  if(upper_limit_mzrange > 0){
    #print('scaled')
    p=p+xlim(NA,upper_limit_mzrange)
  }
# make outputs
  if(print_plot){
    print(p)
  }
  if(return_plot){
    return(p)
  }else if(return_dataframe){
    return(df)
  }else{
    return(NULL)
  }
}
```
[back to top](#table-of-contents)

### msp-Viewer
This is a Shiny app that allows to conveniently open msp files, convert it to table and view the contained mass spectra. So it doesn't need NISTMS or similar tools.
```R
library(shiny)
library(shinythemes)
library(tidyverse)
library(plotly)

# Define UI for application that draws a histogram
ui = fluidPage(
    theme = shinythemes::shinytheme("spacelab"),

 titlePanel("Mass-Spectra - Viewer"),
    
 sidebarLayout(
   sidebarPanel(
       h4("This app reads mass spectra files (.msp) and allows an easy overview about the content."),
       h4("The compound and peak list can be exported and the mass spectra can checked."),
       br(),
       fileInput("InputFile", "Input file",
                 placeholder = "need to be in .msp-format"),
       br(),
       checkboxInput("SplitEntities", "Split Name into different fields? 
                     Switch off when using different databse-style.", value = TRUE),
       textInput("sep_entries", "Separator between entries in peak list",
                 value = ";"),
       textInput("sep_mass_intens", "Separator between m/z and intensity in peak list",
                 value = " "),
       h4(textOutput("Stats")),
       downloadButton("DL_Table",
                      "Download Table as .csv")

        ),   # Sidebarpanel
        
  mainPanel(
        tabsetPanel(
          tabPanel("Table",
                   br(),
                   DT::dataTableOutput("Table")  ), 
          tabPanel("Spectra",
                   br(),
             fluidRow(
             column(width = 7,
                    uiOutput("CompoundSelector")) #,
    #         column(width = 5,
    #             selectInput("Replicate", "Select one of the replicates to show",
    #                         choices = c("1","2","3","4","5")))
             ),
             fluidRow(
             column(width = 4,
                    sliderInput("topMasses", "Show top N masses:", 
                                min = 0, max = 30, value = 10, step =1),
                    h4("Comments:"),
                    verbatimTextOutput("Comments")),
             column(width = 8,
                    plotlyOutput("Spectra"))
             )) # tab sepctra

            
        )) # Main Panel
    ) # sidebarlayout
) # UI

# Define server logic required to draw a histogram
server = function(input, output) {
    
    # Imort data
    InputFile = reactive({ 
        infile <- input$InputFile
        if (is.null(infile)) {  
            return(NULL) 
        } else {
            return(read_msp(input$InputFile$datapath))        
        } 
     })
    
    # show base stats for file
    output$Stats = renderText({
        req(InputFile())
        paste0("The uploaded dataset contains ", length(unique(InputFile()$Name)),
               " unique compounds with ", length(InputFile()$Peaks),
               " spectra.")
    })
    
    # processed table - split Name into identifier
    # this works only when the Kempa-Ident naming rules are fulfilled
    # this can also be a good check for consistency!
    FinalTable = reactive({
      req(InputFile())
      
      if (input$SplitEntities == TRUE) {
        FinalTable = InputFile() %>% 
            mutate(Name = str_remove(Name, "RI:"),
                   Name = str_remove(Name, regex("Ident:", ignore_case = TRUE))) %>% 
            separate(Name, into = c("Compound", "Derivative", "Product", "RI", "Ident"),
                     sep = "_")
      } else FinalTable = InputFile()
    })

    # results table
    output$Table = DT::renderDataTable(FinalTable())
    
    # download results
    output$DL_Table = downloadHandler(
        filename = "Table.csv",
        content = function(file) {
            write_csv(FinalTable() , file) },
        contentType = "file/csv" )
    
    # compound list for the selection
    CompoundList = reactive({
        req(InputFile())
        sort(InputFile()$Name)
        })
    
    output$CompoundSelector = renderUI({
        req(InputFile())
        selectInput("Compound", "Select the compound to show",
                    CompoundList(), multiple = FALSE,
                    width = "100%")
    })

    # filtered table
    Results = reactive({
        req(InputFile())
        req(input$Compound)
    
    filter(InputFile(), Name == input$Compound)
    })
    
    # show comments for selected compound
    output$Comments = renderText({
        req(Results())
            pull(Results(), Comments)
    })
    
    output$Spectra = renderPlotly({
        req(Results())
        
      spectra_sep = paste0(input$sep_mass_intens, "+")  
      
      #sep_mass_intens
        spectra_df = Results() %>% 
            select(Peaks) %>% 
            separate_rows(Peaks, sep = input$sep_entries) %>%
            #mutate(Data = str_trim(Data, side = "both")) %>% 
            separate(Peaks, into = c("Mass", "Intensity"),
                     sep = spectra_sep,
                     convert = TRUE) %>% 
            filter(!is.na(Intensity)) %>% 
            arrange(desc(Intensity)) %>% 
            mutate(Rank = rank(desc(Intensity), ties.method = "first" ))
        
        ggplotly(ggplot(spectra_df, aes(x= Mass, y = Intensity)) +
            geom_col(width = 0.8) + 
            theme_classic(base_size = 16) + 
            geom_text(data = filter(spectra_df, Rank <= input$topMasses), 
                      aes(label = Mass, y = Intensity + 30), 
                      check_overlap = TRUE) + 
            labs(x = "m/z",
                 title = paste0("\n", input$Compound)) +
          theme(plot.title = element_text(hjust = 0.5,size = 12))
        )
        })
    
} # Server 

# Define the functions:
read_msp = function(inputfile) {

    ## Read in the msp
    input = read_file(inputfile) %>%
        str_replace_all("\r*\n", "\n") %>% 
        str_remove("\n\n$")   # remove empty line at the end
    
    ## transform the file for better readability:
    # split the different entries, separated by a blank line
    input2 = str_split(input, "\n\n", simplify = FALSE)  # generates list
    
    input3 = unlist(input2)                                  # back to character vector
    
    ## extract the important information
    # extract compound names: everything between "Name:" and the next line break
    names = str_extract(input3, regex("(?<=^(Name|NAME|name): ).*(?=\n)", # somehow ignore case doesnt work?
                        ignore_case = TRUE))  
    
    # the next steps aren't very elegant but work
    # remove the names from the input (as this nicely ends with a line break)
    # line breaks mess up with the extraction of the other features, 
    # so replace line-breaks by "___" to restore them later 
    
    input_wo_names = str_remove(input3, regex("^(Name|NAME|name): .*\n", ignore_case = TRUE)) %>%  
        str_replace_all("\n","___")
    
    # extract comments: everything from after the name and (including "Num Peaks: NNN)
    # then restore line breaks
    comments = str_replace_all(str_extract(input_wo_names,
                                           regex(
                                           ".*Num Peaks: *\\d{1,3}", ignore_case = TRUE)), 
                               "___", "\n")
    
    #extract peaks: everything after "Num Peaks: NNN" to the end
    peaks = str_extract(input_wo_names,
                        regex("(?<=Num Peaks: \\d{1,3}___).*" , ignore_case = TRUE))
    # clean up, remove unneeded symbols
    peaks = str_replace_all(peaks,
                            c("___" = "",     # old line breaks
                              "; *" = ";",    # empty space after semicolon
                              "^ *" = "",     # empty space in front of first number
                              ";$" = ""))     # ending semicolon in last row
    
    # create a dataframe by combining the 3 vectors:
    msp_as_df = tibble("Name" = names,
                       "Comments" = comments,
                       "Peaks" = peaks) %>% 
      #mutate(RI=str_replace_all(Name,".*RI:(\\d+)_.*","\\1")) %>% 
      select(Name,Comments,Peaks)
    
    # check the results:
    if (is.data.frame(msp_as_df) & all(colnames(msp_as_df) == c("Name", "Comments", "Peaks")) == TRUE) {
        print(paste0("Done! ", length(unique(msp_as_df$Name)),
                     " unique compounds with ", length(msp_as_df$Peaks),
                     " spectra imported"))
        
        return(msp_as_df)
    } else {
        print("Error, conversion failed!")
        return(msp_as_df) }
    
}  # Function read_msp


# Run the application 
shinyApp(ui = ui, server = server)
```
[back to top](#table-of-contents)

# A collection of R snippets

## Table of Contents

* [Metabolomics](#metabolomics)
   * [Data processing](#data-processing)
      * [Split long metabolite names](#Split-long-metabolite-names)
   * [Handling .msp-files](#handling-.msp-files)
      * [Read msp-files](#read-msp-files)
      * [Plot Spectra from msp-data style](#plot-msp)
      * [msp-Viewer](#msp-viewer)

## Metabolomics

### Data processing

This section reports some potential tools to process the data

#### Split long metabolite names

```R
  ###   This script illustrate different approaches to focus on the relevant metabolite information
  ###   obtained from metaolomics experiments.
  ###   Typically (e.g. for plotting) you don't want to show the full name of the database entry
  ###   but just the metabolite name, maybe with the derivatisation and mainproduct/ byproduct information
  ###    

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

  ## Now we want to simplify the name, for this we have multiple options available: 

  ## OPTION 1 - Just take the name as entry before the second underscore  - Thanks to Tobi
  data_sep1A = mutate(data, Name = sub("(.*?_.*?)_(.*)","\\1", Compound))

  # doing the same using str_extract() in the stringr library (tidyverse)
  library(stringr)
  data_sep1B = mutate(data, Name = str_extract(Compound, ".*?_.*?(?=_.*)"))

  data_sep1A == data_sep1B

  ## OPTION 2 - Just split into all entities along the underscore
  # this might lead to problems, when there are missing or too much underscores, see for example Leucine
   data_sep2 = separate(data, Compound, into = c("Name", "Derivative",
                                                "Product", "RI", "Ident"), sep = "_") 

  # 2B - splitted data can then merged together again
  data_united = unite(data_sep2, Metabolite, c("Name", "Derivative", "Product"), sep = "_") 


  ## OPTION 3. Split at the second underscore 
   data_sep3 = 


  ## OPTION 4 - FAIL_SAFE APPROACH: Split before the RI:
  # Clearly this only works when "_RI:" is there   
   data_sep4 = separate(data, Compound, into = c("Name", "Details"), sep = "_(?=RI:)") 
```

### Handling .msp-files

This section reports some ideas to work with mass spectra (msp) files 

#### Read msp-files

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

#### Plot Spectra from msp-data style

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

plotspectra=function(stringfromMSP, mytitle="substance", print_plot=T, return_plot=F, return_dataframe=F, normalize_Intensity=TRUE,
                     upper_limit_mzrange=0){
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



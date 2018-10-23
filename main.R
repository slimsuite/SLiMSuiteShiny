################# ::: APP INFO ::: ######################
info = list(
  apptitle = "SLiMSuite REST Parser",
  version = "0.1.0",
  lastedit = "29/09/2017",
  author = "Richard J. Edwards",
  contact = "richard.edwards@unsw.edu.au",
  description = "Generic SLiMSuite Shiny app for parsing [SLiMSuite REST](http://rest.slimsuite.unsw.edu.au/) output."
)
#i# The main.R script load libraries and contains the initial parameter settings and functions.
#i# This script is used by the `ui.R` and `server.R` files.

################# ::: HISTORY ::: ######################
# V0.1.0 - Development version based on Template V0.1.0.
# V0.2.0 - Added extra HTML output to add reactive Results tab info.

################## ::: LICENSE ::: #####################
# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# SLiMEnrich program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License (gpl.md)
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

################# ::: KNOWN BUGS ::: ######################
#!# List bugs here.

############### ::: TO DO ::: ##################
# https://stackoverflow.com/questions/19470426/r-shiny-add-tabpanel-to-tabsetpanel-dynamically-with-the-use-of-renderui
# https://stackoverflow.com/questions/35020810/dynamically-creating-tabs-with-plots-in-shiny-without-re-creating-existing-tabs/

############### ::: GENERAL SETUP ::: ##################
#i# This used to be the setup.R code...

############### ::: REQUIRED LIBRARIES ::: ##################
# Check whether packages of interest are installed
is_installed = function(mypkg) is.element(mypkg, installed.packages()[,1]) 
# Install library if not already installed
# Run a for-loop of all the package names listed below in the function call
# with the list of packages: load_or_install(c("pkg1", "pkg2",..., "pkgn"))
load_or_install = function(package_names) 
{ 
  for(package_name in package_names) 
  { 
    if(!is_installed(package_name)) 
    { 
      #install.packages(package_name,repos="http://lib.stat.cmu.edu/R/CRAN") 
      install.packages(package_name)
    } 
    library(package_name,character.only=TRUE,quietly=TRUE,verbose=FALSE) 
  } 
}
# Load or install required libraries
load_or_install(c("shiny", "httr", "DT", "markdown","plyr","tools","shinyjs","igraph","visNetwork"))

############### ::: SET DEFAULTS ::: ##################
settings = list(
  #i# URL of REST servers
  resturl = "http://rest.slimsuite.unsw.edu.au/",
  #i# URL of Dev REST servers
  devurl = "http://restdev.slimsuite.unsw.edu.au/",
  #i# Standard REST server &rest keys
  stdkeys = c("status", "version", "ini", "log", "warnings", "errors", "outfmt", "help", "jobid", "intro", "prog"),
  #i# Subset of REST server keys that are not returned by restKeys()
  restkeys = c("errors", "outfmt", "help", "intro", "prog", "compare","self-compare"),
  #i# Define REST output formats
  outfmts = c("none","text","csv","tsv","log"),
  csv = c("main","occ"),
  tsv = c("compare","self-compare"),
  #i# Whether to run in debugging mode
  debug = TRUE,
  #i# Default program
  prog = "zen",
  #i# Default JobID
  jobid = "17092800011" # "17092300003"   # Test with 17092500002 & 17092800011
)

############### ::: SETUP DATA ::: ##################
#i# This function is called by server.R in a reactiveValues() call.
#i# This way, it will be updated whenever the settings that affect it are changed.
setupData = function(){
  emptydb = list()
  for(rkey in settings$stdkeys){
    if(rkey == "log"){
      emptydb[[rkey]] = data.frame()
    }else{
      emptydb[[rkey]] = ""
    }
  }
  emptydb$intro = "1. Use \"Upload Data\" to get a JobID.\n2. If already have a JobID, use \"Retrieve Job\"."
  # Return data list
  return(emptydb)
}

############### ::: REST FUNCTIONS ::: ##################
### Check whether a jobID looks legit and return True or False
isJobID <- function(jobid){
  if(nchar(jobid) != 11){ return(FALSE) }
  if(is.na(as.numeric(jobid))){ return(FALSE) }
  if(is.na(as.integer(substr(jobid,1,6)))){ return(FALSE) }
  if(is.na(as.integer(substr(jobid,7,11)))){ return(FALSE) }
  return(TRUE)
}

### Check whether a file exists and return True or False
isFile <- function(FASTA){
  if(length(FASTA)==0){return(FALSE)}
  #if(file_ext(FASTA)=="fasta"){return(TRUE)}
  return(TRUE)
}
### Check whether the sequences in the file is leagal
isSequence<-function(seq){
  if(substr(seq,1,1) == '>'){return(TRUE)}
  else{return(FALSE)}
}
### If the sequences are leagal, modify them to meet the requirement to put in url
modifySequence <- function(seq) {
  seq_list = strsplit(seq, "\\n")[[1]]
  modified = ""
  i = 1
  n = 1
  while (n <= length(seq_list)) {
    s = seq_list[n]
    if (s != "") {
      prefix = strsplit(s, " ")[[1]][1]
      # if there is not uniprot ID for this sequence, use the original tag
      if (prefix == "") {
        prefix = s
      }
      seq = paste0(prefix, ":", seq_list[n+1])
    }
    if (modified != "") {
      modified = paste0(modified, ",", seq)
    } else {
      modified = seq
    }
    n = n + 2
  }
  return (modified)
}

### Check whether Job has run
checkJob <- function(jobid,password=""){
  checkurl = paste0(settings$resturl,"check&jobid=",jobid,"&password=",password)
  jobcheck = readLines(checkurl,warn=FALSE)[1]
  if(jobcheck == "Finished"){ return(TRUE) }
  else{ return(jobcheck) }  # Missing/Queued/Running
  #i# Password mismatch: "ERROR: JobID password mismatch! Contact server admin if this is your job and you have forgotten the password used."
} 
### Function for returning the REST keys
getRestKeys <- function(jobid,password=""){
  joburl = paste0(settings$resturl,"retrieve&jobid=",jobid,"&rest=restkeys")
  return(readLines(joburl,warn=FALSE))
}
### Return the REST output format
getRestFormat <- function(rest,pure=TRUE){
  outfmt = "text"
  if(rest %in% c("None","none","")){ return("none") }
  if(rest %in% settings$csv){ outfmt = "csv" }
  if(rest %in% settings$tsv){ outfmt = "tsv" }
  if(rest == "log"){ outfmt = "log" }
  if(pure == FALSE & outfmt %in% c("csv","tsv","log")){ outfmt = "table" }
  return(outfmt)
}
### Return an R object with REST output
getRestOutput <- function(jobid,rest,outfmt="",password=""){
  joburl = paste0(settings$resturl,"retrieve&jobid=",jobid,"&rest=",rest)
  if(outfmt == ""){
    outfmt = getRestFormat(rest)
  }
  if(outfmt == "text"){ return(readLines(joburl,warn=FALSE)) }
  if(outfmt == "csv"){ 
    # add CompariMotifs link in main table
    if (rest == "main"){
      motifs <- read.delim(joburl,header=TRUE,sep=",",stringsAsFactors=FALSE)
      links <- c()
      for (i in 1:(length(motifs$Pattern))){
        link <- paste0("http://rest.slimsuite.unsw.edu.au/comparimotif&motifs=",motifs$Pattern[i],"&searchdb=elm")
        text <- c("Show CompariMotif Information")
        links <- c(links,HTML("<a href=",link,">",text,"</a>"))
      }
      motifs$Link <- links
      return(motifs)
    }
    # add ProViz link in occ table
    if (rest == "occ"){
      motifs <- read.delim(joburl,header=TRUE,sep=",",stringsAsFactors=FALSE)
      links <- c()
      for (i in 1:(length(motifs$Seq))){
        uniprotID <- strsplit(motifs$Seq[i], "__")[[1]][2]
        # check if what we get is in uniprot ID format: P22363
        isUniprotID = grep("\\w\\d+", uniprotID)
        # only generate url for having Uniprot ID
        if (length(isUniprotID) > 0 && isUniprotID == 1) {
          start <- motifs$Start_Pos[i]
          end <- motifs$End_Pos[i]
          showStart = start - 40
          showEnd = end + 40
          pattern_match <- motifs$Match[i]
          link <- paste0("http://proviz.ucd.ie/proviz.php?uniprot_acc=",uniprotID,"&ali_start=",showStart,"&ali_end=",showEnd,
                         "&tracks=peptides,SLiMFinder,AAA,",pattern_match,",",start,",",end)
          text <- c("Show ProViz Information")
          links <- c(links,HTML("<a href=",link,">",text,"</a>"))
        } else {
          links <- c(links, HTML("Didn't find Uniprot ID in given file. Cannot generate ProViz information without Uniprot ID."))
        }
      }
      motifs$ProViz <- links
      return(motifs)
    }
    return(read.delim(joburl,header=TRUE,sep=",",stringsAsFactors=FALSE))}
  if(outfmt == "tsv"){ 
    #joburl = paste0(settings$resturl,"retrieve&jobid=","18092800018","&rest=",rest)
    return(read.delim(joburl,header=TRUE,sep="\t",stringsAsFactors=FALSE)) }
  if(outfmt == "log"){ 
    logdata = read.delim(joburl,header=FALSE,sep="\t")
    # Modify Log data
    colnames(logdata) = c("Type","Time","Details")
    logdata$Line = as.integer(rownames(logdata))
    logdata = logdata[,c(4,1,2,3)]
    return(logdata) 
  }
  #if (outfmt == 'plot'){
  #  if (rest== 'cloud'){return(read.delim(joburl,header=TRUE,sep=",",stringsAsFactors=FALSE))}}
}
### Check the masking options
# http://rest.slimsuite.unsw.edu.au/docs&page=module:qslimfinder
CheckMask <- function(dismask,consmask,ftmask,imask){
  if (dismask == TRUE){joburl_d = paste0(settings$resturl,"slimfinder&dismask=","T")}
  else{joburl_d = paste0(settings$resturl,"slimfinder&dismask=","F")}
  
  if (consmask == TRUE){joburl_c = paste0(joburl_d,"&consmask=","T")}
  else{joburl_c = paste0(joburl_d,"&consmask=","F")}
  
  if (ftmask == "EM"){joburl_f = paste0(joburl_c,"&ftmask=","EM")}
  else if (ftmask=="DOMAIN"){joburl_f = paste0(joburl_c,"&ftmask=","DOMAIN")}
  else if (ftmask=="TRANSMEM") {joburl_f = paste0(joburl_c,"&ftmask=","TRANSMEM")}
  else{joburl_f =  paste0(joburl_c,"&ftmask=","")}
  
  if (imask == "inclusively"){joburl_i = paste0(joburl_f,"&imask=","inclusively")}
  else {joburl_i = paste0(joburl_f,"&imask=","")}
  
  return(joburl_i)
}

### Return a JobID(Input: sequences) with REST output
getSequences <- function(se,dismask,consmask,ftmask,imask){
  maskoptions <- CheckMask(dismask,consmask,ftmask,imask)
  joburl = paste0(maskoptions,"&seqin=",se)
  result <- readLines(joburl,warn=FALSE)
  return(substr(result[96], 20,30)) 
}

### Return a JobID(Input: UniProtID) with REST output
getUniprotID <- function(id,dismask,consmask,ftmask,imask){
  maskoptions <- CheckMask(dismask,consmask,ftmask,imask)
  joburl = paste0(maskoptions,"&uniprotid=",id)
  result <- readLines(joburl,warn=FALSE)
  return(substr(result[96], 20,30)) 
}

### Return a CompareMotif ID according to jobID
getCompareID <- function(id){
  url_1 = paste0(settings$resturl,"comparimotif&motifs=rest:",id,":main&searchdb=elm")
  #  check the status of the job
  status_1 = readLines(url_1,warn=FALSE)
  while((! status_1 %in% c("Finished","Failed","refresh"))&&(status_1 %in%c("refresh"))){
    writeLines(paste0(status_1," - Sleep!"))
    Sys.sleep(5)   # We will pause for a five seconds to give the job a chance
    status_1 = readLines(url_1,warn=FALSE)
  }
  result_1 <- readLines(url_1,warn=FALSE)
  return(substr(result_1[96], 22,32))
}

### Return a Self CompareMotif ID according to jobID
getSelfCompareID <- function(id){
  url_2 = paste0(settings$resturl,"comparimotif&motifs=rest:",id,":main&searchdb=rest:",id,":main")
  #  check the status of the job
  status_2 = readLines(url_2,warn=FALSE)
  while((! status_2 %in% c("Finished","Failed","refresh"))&&(status_2 %in%c("refresh"))){
    writeLines(paste0(status_2," - Sleep!"))
    Sys.sleep(5)   # We will pause for a five seconds to give the job a chance
    status_2 = readLines(url_2,warn=FALSE)
  }
  result <- readLines(url_2,warn=FALSE)
  return(substr(result[96], 22,32))
}

############### ::: PLOT FUNCTIONS ::: ##################
### Return protein nodes for graph
getProteinNodes <- function(jobid){
  joburl = paste0(settings$resturl,"retrieve&jobid=",jobid,"&password=","","&rest=occ")
  x <- read.delim(joburl,header=TRUE,sep=",",stringsAsFactors=FALSE)
  a <- c()
  # find protein nodes
  for (se in x$Seq){
    register <- TRUE
    for (s in a){
      if (identical(se,s)){
        register<- FALSE
      }}
    if (register){
      a<-c(a,se)
    }
  }
  return(a)
}
### Return motif nodes for graph
getMotifNodes <- function(jobid){
  joburl = paste0(settings$resturl,"retrieve&jobid=",jobid,"&password=","","&rest=occ")
  x <- read.delim(joburl,header=TRUE,sep=",",stringsAsFactors=FALSE)
  a <- c()
  # find motif nodes
  for (pa in x$Pattern){
    register <- TRUE
    for (p in a){
      if (identical(pa,p)){
        register<- FALSE
      }}
    if (register){
      a<-c(a,pa)
    }
  }
  return(a)
}

### Return the edges for graph
getEdges<- function(jobid,nodes){
  node1<-c()
  node2<-c()
  # protein&motif relation
  joburl = paste0(settings$resturl,"retrieve&jobid=",jobid,"&password=","","&rest=occ")
  x <- read.delim(joburl,header=TRUE,sep=",",stringsAsFactors=FALSE)
  for (i in 1:length(x$Pattern)){
    node1<-c(node1,match(x$Seq[i],nodes))
    node2<-c(node2,match(x$Pattern[i],nodes))
  }
  # motif&motif relation
  joburl = paste0(settings$resturl,"retrieve&jobid=",jobid,"&password=","","&rest=main")
  y <- read.delim(joburl,header=TRUE,sep=",",stringsAsFactors=FALSE)
  for (i in 1:length(y$Cloud)){
    for (j in 1:length(y$Cloud)){
      if (y$Cloud[i]==y$Cloud[j]){
        if (i != j){
          node1<-c(node1,match(y$Pattern[i],nodes))
          node2<-c(node2,match(y$Pattern[j],nodes))
        }
      }
    }
  }
  # protein&protein relation
  joburl = paste0(settings$resturl,"retrieve&jobid=",jobid,"&password=","","&rest=upc")
  z <- read.delim(joburl,header=TRUE,sep=",",stringsAsFactors=FALSE)
  for (i in 2:length(z[,1])){
    proteins <- strsplit(z[,1][i], '\t')[[1]][4]
    protein<- strsplit(proteins," ")
    for(m in 1:length(protein[[1]])){
      for(n in 1:length(protein[[1]])){
        if (m!=n){ 
          node1<-c(node1,match(protein[[1]][m],nodes))
          node2<-c(node2,match(protein[[1]][n],nodes))}
      }
    }
  }
  M <- cbind(node1,node2)
  return(M)
}


############### ::: UPDATE DATA ::: ##################
#!# This is the old function that needs updating with above functions
#i# This function is called by server.R in a reactiveValues() call.
#i# This way, it will be updated whenever the settings that affect it are changed.
setData = function(jobid,prog="retrieve",password="",extra=c(),formats=c()){
  #># jobid = REST server job ID
  #># prog = REST program
  #># password = REST password
  #># extra = Additional 
  restbase = paste("http://rest.slimsuite.unsw.edu.au/",prog,"&jobid=",jobid,"&password=",password,sep="")
  
  # Template SLiMSuite shiny app using example Zen run:
  # http://rest.slimsuite.unsw.edu.au/retrieve&jobid=17092300003&rest=format&password=None&refresh=2
  #i# Standard &rest=X outputs:
  standard = c("status", "version", "ini", "log", "warnings", "errors", "outfmt", "help", "jobid", "intro", "prog")
  #i# Zen-specific &rest=X outputs
  noheadtdt = c("wisdoms")
  headtdt = c()
  noheadcsv = c()
  headcsv = c()
  
  #<# Will return a list of the different REST outputs
  #i# See: https://cran.r-project.org/web/packages/httr/vignettes/quickstart.html
  fulldata = read.delim(paste(restbase,"&rest=full",sep=""),header=FALSE)
  rdata = list(full=fulldata)
  rdata$full = as.character(rdata$full[1:3,1])
  for(rest in c(standard,noheadtdt)){
    #resturl = paste(restbase,"&rest=",rest,sep="")
    rdata[[rest]] <- read.delim(paste(restbase,"&rest=",rest,sep=""),header=FALSE,sep="\t")
  }
  for(rest in c(headtdt)){ rdata[[rest]] <- read.delim(paste(restbase,"&rest=",rest,sep=""),header=TRUE,sep="\t") }
  for(rest in c(noheadcsv)){ rdata[[rest]] <- read.delim(paste(restbase,"&rest=",rest,sep=""),header=FALSE,sep=",") }
  for(rest in c(headcsv)){ rdata[[rest]] <- read.delim(paste(restbase,"&rest=",rest,sep=""),header=TRUE,sep=",") }
  # Modify Log data
  colnames(rdata$log) = c("Type","Time","Details")
  rdata$log$Line = as.integer(rownames(rdata$log))
  rdata$log = rdata$log[,c(4,1,2,3)]
  # Return data list
  return(rdata)
}
#cat(as.character(pep[1:10,1]))



############### ::: SHINY CODE INFO ::: ##################
#i# This section contains some information comments that can be deleted in actual Apps.

# Template SLiMSuite shiny app using example Zen run:
# http://rest.slimsuite.unsw.edu.au/retrieve&jobid=17092300003&rest=format&password=None&refresh=2
#i# Standard &rest=X outputs:
# - status, version, ini, log, warnings
#i# Zen-specific &rest=X outputs
# - wisdoms

# # status: ./17092300003.status
# JobID 17092300003 (zen) Finished.
# IP:58.178.74.181
# No queue.
# zen&wisdoms=10
# Run Started: 2017-09-23 13:37:10; PID=5456
# Run finished: 2017-09-23 13:37:11

#i# Output options
#i# Render text inside <pre> html block:
# verbatimTextOutput(outputId, placeholder = FALSE) 

# htmlOutput (uiOutput) = Create an HTML output element
# plotOutput (imageOutput) = Create an plot or image output element
# outputOptions = Set options for an output object.
# tableOutput (dataTableOutput) = Create a table output element
# textOutput = Create a text output element
# verbatimTextOutput = Create a verbatim text output element
# downloadButton (downloadLink) = Create a download button or link
# Progress = Reporting progress (object-oriented API)
# withProgress (setProgress, incProgress) = Reporting progress (functional API)
# modalDialog = Create a modal dialog UI
# urlModal = Generate a modal dialog that displays a URL
# showModal (removeModal) = Show or remove a modal dialog
# showNotification (removeNotification) = Show or remove a notification

# renderPlot = Plot Output
# renderText = Text Output
# renderPrint = Printable Output
# renderDataTable = Table output with the JavaScript library DataTables
# renderImage = Image file output
# renderTable = Table Output
# renderUI = UI Output
# downloadHandler = File Downloads

#i# UI Inputs:
# actionButton (actionLink) = Action button/link
# checkboxGroupInput = Checkbox Group Input Control
# checkboxInput = Checkbox Input Control
# dateInput = Create date input
# dateRangeInput = Create date range input
# fileInput = File Upload Control
# numericInput = Create a numeric input control
# radioButtons = Create radio buttons
# selectInput (selectizeInput) = Create a select list input control
# sliderInput (animationOptions) = Slider Input Widget
# submitButton = Create a submit button
# textInput = Create a text input control
# textAreaInput = Create a textarea input control
# passwordInput = Create a password input control
# modalButton = Create a button for a modal dialog
# updateActionButton = Change the label or icon of an action button on the client
# updateCheckboxGroupInput = Change the value of a checkbox group input on the client
# updateCheckboxInput = Change the value of a checkbox input on the client
# updateDateInput = Change the value of a date input on the client
# updateDateRangeInput = Change the start and end values of a date range input on the client
# updateNumericInput = Change the value of a number input on the client
# updateRadioButtons = Change the value of a radio input on the client
# updateSelectInput (updateSelectizeInput) = Change the value of a select input on the client
# updateSliderInput = Change the value of a slider input on the client
# updateTabsetPanel (updateNavbarPage, updateNavlistPanel) = Change the selected tab on the client
# insertTab (prependTab, appendTab, removeTab) = Dynamically insert/remove a tabPanel
# showTab (hideTab) = Dynamically hide/show a tabPanel
# updateTextInput = Change the value of a text input on the client
# updateTextAreaInput = Change the value of a textarea input on the client
# updateQueryString = Update URL in browser's location bar
# getQueryString (getUrlHash) = Get the query string / hash component from the URL

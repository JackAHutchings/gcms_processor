# Version 0!

rm(list=ls())

library(dplyr); options(dplyr.summarise.inform = FALSE) # Data processing
library(ggplot2) # Data visualization
library(cowplot); theme_set(theme_classic()) # Additional ggplot functionality
library(lemon) # Additional ggplot functionality
library(tidyr) # Data processing
# library(zoo) # needed for na.locf function
# library(lemon) # needed for facet_rep_wrap/grid functions
library(writexl) # export the finished data
library(readxl) # import the template file
library(shinydashboard) # shiny
library(httr) # Needed for fetching URLs
library(DT) # tabular data is rendered via DT
library(shiny) # shiny
library(shinyFiles) # shiny file functionality
library(plotly) # interactive plots
library(stringr) # data wrangling
# library(htmlwidgets)
# library(shinyjs)
# library(shinyWidgets)
library(admisc) # data wrangling, specifically the numdec function

{
  # install.packages("remotes")
  # remotes::install_github("https://github.com/ethanbass/chromConverter/")
}

library(chromConverter) # read and convert MS data

# Peak finding/fitting functions lifted from Ethan Bass' chromatographR package. (https://ethanbass.github.io/chromatographR/)
{
  find_peaks <- function(y, smooth_type = c("gaussian", "box", "savgol",
                                            "mva", "tmva", "none"),
                         smooth_window = .001, slope_thresh = 0, amp_thresh = 0,
                         bounds = TRUE){
    if (!is.vector(y)){
      stop("Please provide a vector to argument `y` to proceed.")
    }
    smooth_type <- match.arg(smooth_type, c("gaussian", "box", "savgol", "mva",
                                            "tmva", "none"))
    if (smooth_window < 1){
      smooth_window <- max(length(y) * smooth_window, 1)
    }
    
    savgol <- function(T, fl, forder = 4, dorder = 0) {
      stopifnot(is.numeric(T), is.numeric(fl))
      if (fl <= 1 || fl %% 2 == 0)
        stop("Argument 'fl' must be an odd integer greater than 1.")
      n <- length(T)
      
      # -- calculate filter coefficients --
      fc <- (fl-1)/2                          # index: window left and right
      X <- outer(-fc:fc, 0:forder, FUN = "^")   # polynomial terms and coeffs
      Y <- pinv(X);                           # pseudoinverse
      
      # -- filter via convolution and take care of the end points --
      T2 <- convolve(T, rev(Y[(dorder + 1),]), type = "o")   # convolve(...)
      T2 <- T2[(fc+1):(length(T2)-fc)]
      
      Tsg <- (-1)^dorder * T2
      return( Tsg )
    }
    
    # compute derivative (with or without smoothing)
    if (smooth_type == "savgol"){
      if ((smooth_window %% 2) == 0){
        smooth_window <- smooth_window + 1
      }
      d <- savgol(diff(y), fl = smooth_window)
    } else if (smooth_type == "mva"){
      d <- caTools::runmean(diff(y), k = smooth_window)
    } else if (smooth_type == 'gaussian'){
      d <- diff(ksmooth(seq_along(y), y, kernel = "normal",
                        bandwidth = smooth_window)$y)
    }  else if (smooth_type == "box"){
      d <- diff(ksmooth(seq_along(y), y, kernel = "box",
                        bandwidth = smooth_window)$y)
    } else if (smooth_type == "tmva"){
      d <- runmean(runmean(diff(y), k = smooth_window), k = smooth_window)
    } else{
      d <- diff(y)
    }
    
    # detect zero-crossing of first derivative (peak apex)
    p1 <- which(sign(d[1:(length(d)-1)]) > sign(d[2:length(d)]))
    
    # detect second derivative exceeding slope threshold
    p2 <- which(abs(diff(d)) > slope_thresh)
    
    # detect y-vals exceeding amplitude threshold
    p3 <- which(y > amp_thresh) 
    p <- intersect(intersect(p1,p2), p3)
    if (bounds){
      p4 <- which(sign(d[1:(length(d)-1)]) < sign(d[2:length(d)]))
      
      # find lower bound
      suppressWarnings(bl <- sapply(p, function(v) max(p4[p4 < v])))
      bl[which(bl == -Inf)] <- 1
      
      # find upper bound
      suppressWarnings(bu <- sapply(p, function(v) min(p4[p4 > v])))
      bu[which(bu == Inf)] <- length(y)
      p <- data.frame(pos = p, lower = bl, upper = bu)
    }
    p
  }
  
  fit_peaks <- function (resp, #y-axis vector of chromatogram
                         pos = NULL,
                         sd.max = 50,
                         fit = c("egh", "gaussian", "raw"),
                         max.iter = 1000,
                         r2_thresh = 0.5,
                         ...){
    
    fit <- match.arg(fit, c("egh", "gaussian", "raw"))
    if (is.null(pos)){
      pos <- find_peaks(resp, ...)
    }
    
    tabnames <- switch(fit,
                       "gaussian" = c("rt", "start", "end", "sd", "FWHM", 
                                      "height", "area", "r-squared"),
                       "egh" = c("rt", "start", "end", "sd", "tau", "FWHM", 
                                 "height", "area", "r.squared"),
                       "raw" = c("rt", "start", "end", "sd", "FWHM", "height",
                                 "area")
    )
    noPeaksMat <- matrix(rep(NA, length(tabnames)),
                         nrow = 1, dimnames = list(NULL, tabnames))
    on.edge <- sapply(pos$pos, function(x){
      x <= 1 || is.na(resp[x + 1]) || is.na(resp[x - 1])
    })
    pos <- pos[!on.edge,]
    if (nrow(pos) == 0) 
      return(noPeaksMat)
    
    {
      fitpk_gaussian <- function(x, pos, max.iter,noise_threshold = .001, ...){
        y <- x
        xloc <- pos[1]
        peak.loc <- seq.int(pos[2], pos[3])
        suppressWarnings(m <- fit_gaussian(x = peak.loc, 
                                           y = y[peak.loc],
                                           start.center = xloc,
                                           start.height = y[xloc],
                                           max.iter = max.iter)
        )
        area <- sum(diff(peak.loc) * mean(c(m$y[-1], tail(m$y,-1)))) # trapezoidal integration
        r.squared <- try(summary(lm(m$y ~ y[peak.loc]))$r.squared, silent = TRUE)
        c("rt" = m$center, "start" = pos[2], "end" = pos[3], 
          "sd" = m$width, "FWHM" = 2.35 * m$width,
          "height" = y[xloc], "area" = area, "r.squared" = r.squared)
      }
      gaussian <- function(x, center = 0, width = 1, height = NULL, floor = 0){
        
        # adapted from Earl F. Glynn;  Stowers Institute for Medical Research, 2007
        twoVar <- 2 * width * width
        sqrt2piVar <- sqrt( pi * twoVar)
        y <- exp( -( x - center)^2 / twoVar) / sqrt2piVar
        
        # by default, the height is such that the curve has unit volume
        if ( ! is.null (height)) {
          scalefactor <- sqrt2piVar
          y <- y * scalefactor * height
        }
        y + floor
      }
      fit_gaussian <- function(x, y, start.center = NULL,start.width = NULL,start.height = NULL,start.floor = NULL,fit.floor = FALSE,max.iter = 1000){
        # estimate starting values
        who.max <- which.max(y)
        if (is.null(start.center)) start.center <- x[who.max]
        if (is.null(start.height)) start.height <- y[who.max]
        if (is.null(start.width)) start.width <- sum( y > (start.height/2)) / 2
        
        # call the Nonlinear Least Squares, either fitting the floor too or not
        controlList <- nls.control(maxiter = max.iter, minFactor = 1/512,
                                   warnOnly = TRUE)
        starts <- list("center" = start.center, "width" = start.width,
                       "height" = start.height)
        if (!fit.floor) {
          nlsAns <- try(minpack.lm::nlsLM(y ~ gaussian(x, center, width, height),
                                          start = starts, control = controlList), silent = TRUE)
        } else{
          if (is.null(start.floor)) start.floor <- quantile(y, seq(0, 1, 0.1))[2]
          starts <- c(starts, "floor" = start.floor)
          nlsAns <- try(minpack.lm::nlsLM( y ~ gaussian(x, center, width, height, floor),
                                           start = starts, control = controlList), silent = TRUE)
        }
        
        # package up the results to pass back
        
        if (inherits(nlsAns, "try-error")){
          yAns <- gaussian(x, start.center, start.width, start.height, start.floor)
          out <- list("center" = start.center, "width" = start.width,
                      "height" = start.height,
                      "y" = yAns, "residual" = y - yAns)
          floorAns <- if (fit.floor) start.floor else 0
        } else {
          coefs <-coef(nlsAns)
          out <- list( "center" = coefs[1], "width" = coefs[2], "height" = coefs[3],
                       "y" = fitted(nlsAns), "residual" = residuals(nlsAns))
          floorAns <- if (fit.floor) coefs[4] else 0
        }
        if (fit.floor) {
          out <- c( out, "floor" = floorAns)
        }
        
        return(out)
      }
      
      fitpk_egh <- function(x, pos, max.iter, noise_threshold = .001){
        y <- x
        xloc <- pos[1]
        peak.loc <- seq.int(pos[2], pos[3])
        suppressWarnings(m <- fit_egh(peak.loc, y[peak.loc], start.center = xloc,
                                      start.height = y[xloc], max.iter = max.iter)
        )
        r.squared <- try(summary(lm(m$y ~ y[peak.loc]))$r.squared, silent = TRUE)
        # trapezoidal integration
        area <- sum(diff(peak.loc) * mean(c(m$y[-1], tail(m$y, -1))))
        c("rt" = m$center, "start" = pos[2], "end" = pos[3], 
          "sd" = m$width, "tau" = m$tau, "FWHM" = 2.35 * m$width,
          "height" = y[xloc], "area" = area, "r.squared" = r.squared)
      }
      egh <- function(x, center, width,  height, tau, floor = 0){
        result <- rep(0, length(x))
        index <- which(2*width^2 + tau*(x-center)>0)
        result[index] <- height*exp(-(x[index] - center)^2/(2*width^2 + tau*(x[index] - center)))
        return(result)
      }
      fit_egh <- function(x1, y1, start.center = NULL, start.width = NULL,
                          start.tau = NULL, start.height = NULL, start.floor = NULL,
                          fit.floor = FALSE, max.iter = 1000) {
        
        # try to find the best egh to fit the given data
        
        # make some rough estimates from the values of Y
        who.max <- which.max(y1)
        if (is.null(start.center)){
          start.center <- x1[who.max]
        }
        if (is.null(start.height)){
          start.height <- y1[who.max]
        }
        if (is.null(start.width)){
          start.width <- sum(y1 > (start.height/2)) / 2
        }
        if (is.null(start.tau)){
          start.tau <- 0
        }
        # call the Nonlinear Least Squares, either fitting the floor too or not
        controlList <- nls.control(maxiter = max.iter, minFactor = 1/512,
                                   warnOnly = TRUE)
        starts <- list("center" = start.center, "width" = start.width, 
                       "height" = start.height, "tau" = start.tau)
        if (!fit.floor){
          nlsAns <- try(minpack.lm::nlsLM(y1 ~ egh(x1, center, width, height, tau),
                                          start = starts, control = controlList), silent = TRUE)
        } else{
          if (is.null( start.floor)) start.floor <- quantile( y1, seq(0, 1, 0.1))[2]
          starts <- c(starts, "floor" = start.floor)
          nlsAns <- try(minpack.lm::nlsLM(y1 ~ egh(x1, center, width, height, tau, floor),
                                          start = starts, control = controlList), silent = TRUE)
        }
        
        # package up the results to pass back
        if (inherits(nlsAns, "try-error")) {
          yAns <- egh(x1, start.center, start.width, start.height,
                      start.tau, start.floor)
          out <- list("center" = start.center, "width" = start.width,
                      "height" = start.height, "tau" = start.tau,
                      "y" = yAns, "residual" = y1 - yAns)
          floorAns <- if (fit.floor) start.floor else 0
        } else {
          coefs <-coef(nlsAns)
          out <- list( "center" = coefs[1], "width" = coefs[2], "height" = coefs[3],
                       "tau" = coefs[4], "y" = fitted(nlsAns),
                       "residual" = residuals(nlsAns))
          floorAns <- if (fit.floor) coefs[5] else 0
        }
        
        if (fit.floor) {
          out <- c( out, "floor" = floorAns)
        }
        return(out)
      }
      fitpk_raw <- function(x, pos, max.iter, noise_threshold = .001){
        y <- x
        xloc <- pos[1]
        peak.loc <- seq.int(pos[2], pos[3])
        
        # perform trapezoidal integration on raw signal
        area <- sum(diff(peak.loc) * mean(c(y[peak.loc][-1], tail(y[peak.loc],-1))))
        c("rt" = pos[1], "start" = pos[2], "end" = pos[3], 
          "sd" = pos[3] - pos[2], "FWHM" = 2.35 * pos[3] - pos[2],
          "height" = y[xloc], "area" = area)
      }
      }
    fitpk <- switch(fit,
                    "gaussian" = fitpk_gaussian,
                    "egh" = fitpk_egh,
                    "raw" = fitpk_raw)
    huhn <- data.frame(t(apply(pos, 1, fitpk, 
                               x = resp,
                               max.iter = max.iter,
                               noise_threshold = noise_threshold)))
    colnames(huhn) <- tabnames
    huhn <- data.frame(sapply(huhn, as.numeric, simplify = FALSE))
    if (!is.null(sd.max)) {
      huhn <- huhn[huhn$sd < sd.max, ]
    }
    huhn$id <- 1:length(huhn$rt)
    if (fit %in% c("gaussian","egh")){
      huhn <- huhn %>% filter(r.squared > r2_thresh)
    }
    x <- try(huhn[huhn$rt > 0,], silent = TRUE)
    if(inherits(x, "try-error")) NA else x
  }
  
}

#testing
{
  gcms_template_file = "C:/Box/Konecky Lab/Data/Agilent GC-MS (Bradley Lab Cyborg)/pared_down_data/gcms_template.xlsx"
}

# UI
{
    ui <- dashboardPage(
        skin = "green",
        dashboardHeader(title = "GC-MS Postprocessor"),
        dashboardSidebar(
            sidebarMenu(
                menuItem("Ingest Data",tabName = "ingest_tab",icon = icon("table")),
                menuItem("Chromatogram Exploration",tabName = "chromatogram_tab",icon = icon("chart-line"))
            )
        ),
        dashboardBody(
            tags$style("* {font-family: Arial;}"),
            tags$h1(tags$style("h1 {font-family: Arial;}")),
            tags$h3(tags$style("h3 {font-family: Arial;}")),
            tabItems(
                tabItem(tabName = "ingest_tab",
                        fluidPage(
                            box(title = h1("Shiny GCMS Quantification App ", strong("(expand for instructions!)"),hr()),
                                solidHeader=T,
                                width=12,
                                collapsible = T,
                                collapsed = T,
                                uiOutput("introduction")),
                            box(title = "Select GC-MS template:",
                                width=3,
                                shinyFilesButton(id="gcms_template",
                                                 title="This must follow the structure of the original template! Do not delete any sheets!",
                                                 multiple=F,
                                                 label = "Browse...")),
                            column(3,
                                infoBoxOutput("gcms_template_status",width=12)),
                            box(title = "Troubleshooting",
                                width= 12,
                                verbatimTextOutput("troubleshooting")),
                            box(title = "Errors",
                                width=12,
                                verbatimTextOutput("errors")),
                            box(title = "Raw GCMS Data",
                                width=12,
                                collapsible = T,
                                collapsed = T,
                                DTOutput("raw_gcms_data"),
                                style="height:500px; overflow-y: scroll;overflow-x: scroll"),
                            box(title = "Sample Info",
                                width=12,
                                collapsible = T,
                                collapsed = T,
                                DTOutput("sample_info"),
                                style="height:500px; overflow-y: scroll;overflow-x: scroll"),
                            box(title = "Standard Info",
                                width=12,
                                collapsible = T,
                                collapsed = T,
                                DTOutput("standard_info"),
                                style="height:500px; overflow-y: scroll;overflow-x: scroll"),
                            box(title = "sequence (direct read from template)",
                                width=12,
                                collapsible = T,
                                collapsed = T,
                                DTOutput("sequence"),
                                style="height:500px; overflow-y: scroll;overflow-x: scroll"),
                            box(title = "comp_list (direct read from template)",
                                width=12,
                                collapsible = T,
                                collapsed = T,
                                DTOutput("comp_list"),
                                style="height:500px; overflow-y: scroll;overflow-x: scroll"),
                            box(title = "parameters (direct read from template)",
                                width=12,
                                collapsible = T,
                                collapsed = T,
                                DTOutput("parameters"),
                                style="height:500px; overflow-y: scroll;overflow-x: scroll"),
                            box(title = "GCMS Line Names",
                                width=12,
                                collapsible = T,
                                collapsed = T,
                                DTOutput("gcms_names"),
                                style="height:500px; overflow-y: scroll;overflow-x: scroll")
                        )
                ),
                tabItem(tabName = "chromatogram_tab",
                        fluidPage(
                          fluidRow(
                            box(title = "Injections to View",
                                width = 4,
                                style='overflow-x: scroll;height:200px;overflow-y: scroll;',
                                collapsible = T,
                                checkboxGroupInput("injection_selection",
                                                   label = NULL,
                                                   choices ="No data loaded")
                            ),
                            box(title="Plot Limits",
                                width = 4,
                                collapsible=T,
                                collapsed=T,
                                numericInput("explorer_min_rt",
                                             label = "RT Start:",
                                             value = NA),
                                numericInput("explorer_max_rt",
                                             label = "RT End:",
                                             value = NA)
                            ),
                            box(title="Ions to Plot",
                                width=4,
                                collapsible=T,
                                collapsed=T,
                                verbatimTextOutput("mz_syntax_description"),
                                textInput("explorer_fs_mz_range",
                                          label = "Full Scan:",
                                          value = NA),
                                textInput("explorer_sim_mz_range",
                                          label = "SIM Scan:",
                                          value = NA)
                            )
                          ),
                          box(title=NULL,
                              width=12,
                              textOutput("explorer_spectra_range"),
                              plotlyOutput("explorer_chromatogram_plot")),
                          box(title=NULL,
                              width=12,
                              plotlyOutput("sample_spectra_plot"))
                        )
                )    
            )
        )
    )
}

server <- function(input, output, session) {
  # Ingest Data Tab
  {
    ingest_function <- function(gcms_template_file) {
      
      # Read in the gcms_template file
      {
        sequence <- read_xlsx(gcms_template_file,sheet="sequence") %>% rename(line=injection)
        comp_list <- read_excel(gcms_template_file,sheet="comp_list")
        parameters <- read_excel(gcms_template_file,sheet="parameters")[,1:2] %>% spread(parameter,value) %>% 
          type.convert(.,as.is=T) %>% # Very neat trick to coerce all columns to the variable type they look like.
          mutate(rfs_spike_volume = ifelse(is.na(rfs_spike_volume),0,rfs_spike_volume))
      }
      
      # Look for ".D" files located in the same directory as the GCMS template file.
      # Convert Hybrid MS Data to excel files containing one sheet each for Full Scan and SIM Scan data. 
      # If exported data already exists, then skip the conversion.
      {
        files <- list.files(path=dirname(gcms_template_file),pattern=".D",full.names=T)
        for(i in 1:length(files)){
          if(!length(list.files(path=files[i],pattern="chrom_data.xlsx")) == 1)
          {
            rawdata <- read_chroms(files[i],format_in=parameters$format_in,format_out="data.frame")
            fullscan <- rawdata[[1]]$MS1 %>% mutate(mz = round(mz,parameters$mz_digits)) %>% group_by(rt,mz) %>% summarize(intensity = sum(intensity))
            simscan <- rawdata[[2]]$MS1
            write_xlsx(list("full"=fullscan,"sim"=simscan),path=paste0(files[i],"/chrom_data.xlsx"))
            rm(rawdata,fullscan,simscan)
          }
        }
        rm(files,i)
      }
      
      # Read back in the converted excel files.
      {
        files <- list.files(dirname(gcms_template_file),pattern="chrom_data.xlsx",recursive=T,full.names=T)
        filenames <- basename(dirname(files))
        
        raw_full <- lapply(files,read_excel,sheet="full") # Read in files
        raw_full <- setNames(raw_full,filenames) # Set injection name
        raw_full <- bind_rows(raw_full,.id="injection") # Convert to a single data.frame
        
        raw_sim <- lapply(files,read_excel,sheet="sim")
        raw_sim <- setNames(raw_sim,filenames)
        raw_sim <- bind_rows(raw_sim,.id="injection")
        
        rawdata <- full_join(raw_full %>% mutate(scan_type = "full"),
                             raw_sim %>% mutate(scan_type = "sim"))
          
        
        rm(files,filenames,raw_full,raw_sim)
      }
      
      # Read and format metadata from the injection folders and template
      {
        # Get entered names for each injection
        {
          files <- list.files(dirname(gcms_template_file),pattern="runstart.txt",recursive=T,full.names=T)
          filenames <- basename(dirname(files))
          
          gcms_names <- lapply(files,readLines) %>% 
            setNames(filenames) %>% 
            bind_rows(.id="injection") %>% 
            gather(injection,runstart) %>% 
            group_by(injection) %>% 
            filter(grepl("dataname",runstart) & !grepl("dataname2",runstart)) %>% 
            mutate(id1 = substr(runstart,gregexpr("= ",runstart)[[1]][1]+2,str_length(runstart))) %>% 
            select(-runstart) %>% 
            mutate(line = as.numeric(substr(injection,gregexpr("_",injection)[[1]][1]+1,gregexpr(".D",injection)[[1]][1]-1)))
          
          rm(files,filenames)
        }
        
        # Grab the standards off the template's sequence sheet, match against chemstation data, and format relevant metadata.
        standard_info <- sequence %>% 
          mutate(id2 = suppressWarnings(as.numeric(id2))) %>% filter(!is.na(id2)) %>% 
          mutate(dummy=T) %>% full_join(parameters %>% mutate(dummy=T),by="dummy") %>% 
          full_join(comp_list %>% mutate(dummy=T),by="dummy",relationship = "many-to-many") %>% 
          left_join(gcms_names,by=c("line","id1")) %>% 
          select(line,injection,id1,id2,comp:mixed_stock_concentration,cal_curve_dilution_factor,cal_curve_volume,rfs_spike_volume,rt_window) %>% 
          mutate(total_volume = cal_curve_volume + rfs_spike_volume,
                 initial_concentration = mixed_stock_concentration * cal_curve_volume / total_volume,
                 initial_concentration = ifelse(is.na(initial_concentration),0,initial_concentration),
                 curve_concentration = initial_concentration * cal_curve_dilution_factor^id2,
                 window_lower = rt-rt_window,
                 window_upper = rt+rt_window
          )
        
        # Do the same thing for samples
        sample_info <- sequence %>% filter(grepl("sample",id2)) %>% 
          mutate(dummy=T) %>% full_join(parameters %>% mutate(dummy=T),by="dummy") %>% 
          left_join(gcms_names,by=c("line","id1")) %>% 
          select(line,injection,id1:sample_mass_grams,rfs_spike_volume) %>% 
          mutate(total_volume = volume + rfs_spike_volume)
      }
      
      # Some simple summary info about our rawdata.
      {
        scan_info <- rawdata %>% 
          group_by(scan_type) %>% 
          mutate(min_mz = min(mz),
                 max_mz = max(mz)) %>% 
          filter(rt < min(rt)+2) %>% 
          group_by(scan_type,min_mz,max_mz) %>% 
          summarize(resolution = numdec(mz)) %>% 
          mutate(mz_string = paste0(min_mz,"-",max_mz))
      }
      
      # Checks for problems
      {
        # check to see if any of our standard or sample entries in the 'sequence' sheet are missing a matching injection.
        {
          injection_sequence_matching <- bind_rows(standard_info %>% select(injection,id1),
                                                   sample_info %>% select(injection,id1)) %>% 
            distinct() %>% 
            mutate(id1 = ifelse(is.na(id1),"",id1)) %>% 
            filter(is.na(injection)) %>% 
            summarize(list = str_c(id1,collapse=" | ")) %>% 
            .$list
          
          
          
          if(injection_sequence_matching > 1){
            injection_sequence_error <- paste0("These id1 entries in the 'sequence' sheet have no matching injection sample name: ",injection_sequence_matching,
                                               "\nSee GCMS Names output box for injection sample names to help find what is wrong.")
          } else {injection_sequence_error = NULL}
          
          
        }
        # check for empty data frames
        {
          empty_list <- data.frame("rawdata" = length(rawdata)==0,
                                   "scan_info" = length(scan_info)==0,
                                   "sample_info" = length(sample_info)==0,
                                   "standard_info" = length(standard_info)==0,
                                   "sequence" = length(sequence)==0,
                                   "comp_list" = length(comp_list)==0,
                                   "parameters" = length(parameters)==0,
                                   "gcms_names" = length(gcms_names)==0
          ) %>% 
            gather(object,empty) %>% 
            filter(empty == T) %>% 
            summarize(list = str_c(object,collapse=" | "))
          
          if(empty_list != ""){
            empty_list_error <- paste0("The following required objects are empty, indicating a problem during import and initial processing: ",empty_list,
                                       "\nCheck that your template and exported GCMS data match examples or contact jackh@wustl.edu")
          } else {empty_list_error = NULL}
          
          errors = c(injection_sequence_error,empty_list_error)
          
        }
        
        final_check = is.null(errors)
        
      }
      
      
      
      if(final_check){
        list("final_check" = final_check,
             "errors" = errors,
             "rawdata" = rawdata,
             "scan_info" = scan_info,
             "sample_info" = sample_info,
             "standard_info" = standard_info,
             "sequence" = sequence,
             "comp_list" = comp_list,
             "parameters" = parameters,
             "gcms_names" = gcms_names)
      } else if(!final_check){
        list("final_check" = errors)
      }
    }
    
    
    
    volumes = c("C"="C:/",getVolumes()())
    
    shinyFileChoose(input,
                    id="gcms_template",
                    roots=volumes,
                    defaultRoot = "C",
                    defaultPath = "/Box/Konecky Lab/Data/Agilent GC-MS (Bradley Lab Cyborg)",
                    session=session,
                    filetypes=c("xlsx"))
    
    observe({
      cat("\ninput$gcms_template value:\n\n")
      print(input$gcms_template)
    })
    
    
    
    output$introduction <- renderUI({
      tagList(
        p("placeholder")
      )
    })
    
    ingest <- reactive({
      if (length(parseFilePaths(volumes,input$gcms_template)$datapath) == 0) {
        return(NULL)
      } else {
        ingest_function(parseFilePaths(volumes,input$gcms_template)$datapath)
      }
    })
    
    output$gcms_template_status <- renderInfoBox({
      if(length(parseFilePaths(volumes,input$gcms_template)$datapath) == 0) {text = "No File Uploaded!"; use_color = "blue"}
      else {text = "File Uploaded."; use_color = "green"}
      infoBox("Status:",text,icon = icon("file-excel"),color = use_color)
    })
    
    output$troubleshooting <- renderPrint(parseFilePaths(volumes,input$gcms_template)$datapath)
    output$errors <-        renderPrint(ingest()$errors)
    output$raw_gcms_data <- renderDT(ingest()$rawdata,options=list('lengthMenu'=JS('[[10,25,50,-1],[10,25,50,"All"]]'),searching=FALSE),class='white-space:nowrap')
    output$sample_info <-   renderDT(ingest()$sample_info,  options=list('lengthMenu'=JS('[[10,25,50,-1],[10,25,50,"All"]]'),searching=FALSE),class='white-space:nowrap')
    output$standard_info <- renderDT(ingest()$standard_info,options=list('lengthMenu'=JS('[[10,25,50,-1],[10,25,50,"All"]]'),searching=FALSE),class='white-space: nowrap')
    output$sequence <-      renderDT(ingest()$sequence,     options=list('lengthMenu'=JS('[[10,25,50,-1],[10,25,50,"All"]]'),searching=FALSE),class='white-space:nowrap')
    output$comp_list <-     renderDT(ingest()$comp_list,    options=list('lengthMenu'=JS('[[10,25,50,-1],[10,25,50,"All"]]'),searching=FALSE),class='white-space:nowrap')
    output$parameters <-    renderDT(ingest()$parameters,   options=list('lengthMenu'=JS('[[10,25,50,-1],[10,25,50,"All"]]'),searching=FALSE),class='white-space:nowrap')
    output$gcms_names <-    renderDT(ingest()$gcms_names,   options=list('lengthMenu'=JS('[[10,25,50,-1],[10,25,50,"All"]]'),searching=FALSE),class='white-space:nowrap')
  }
       
  
  # Chromatogram Explorer
  {
    chromatogram_function <- function(rawdata,
                                      scan_info,
                                      standard_info,
                                      gcms_names,
                                      inj_sel,
                                      min_rt,
                                      max_rt,
                                      full_mz,
                                      sim_mz
                                      ) {
      
      if(!is.null(rawdata)){
        full_mz = ifelse(full_mz=="",paste0(scan_info$min_mz[which(scan_info$scan_type=="full")],"-",
                                            scan_info$max_mz[which(scan_info$scan_type=="full")]),full_mz)
        
        sim_mz = ifelse(sim_mz=="",paste0(scan_info$min_mz[which(scan_info$scan_type=="sim")],"-",
                                          scan_info$max_mz[which(scan_info$scan_type=="sim")]),sim_mz)
        
        mz_selection <- data.frame(scan_type = unique(rawdata$scan_type),
                                   mz_range = c(full_mz,sim_mz)) %>% 
          separate_wider_delim(mz_range,delim=",",too_few = "align_start",names_sep="_") %>% 
          gather(id,range,-scan_type) %>% 
          separate_wider_delim(range,delim="-",too_few="align_start",names_sep="_") %>% 
          gather(range,mz,-c(scan_type,id)) %>% 
          mutate(mz = as.numeric(mz)) %>% 
          spread(range,mz) %>% 
          mutate(range_2 = ifelse(is.na(range_2),range_1,range_2)) %>% 
          full_join(scan_info,by = join_by(scan_type)) %>% 
          select(-c(min_mz,max_mz)) %>% 
          na.omit() %>% 
          rowwise() %>% 
          mutate(mz = list(seq(from=range_1,to=range_2,10^-resolution))) %>% 
          select(scan_type,mz) %>% 
          unnest(mz) %>% 
          distinct()
        
        
        chrom_data <- rawdata %>% left_join(gcms_names) %>% mutate(injection_id1 = paste(injection,id1,sep=" | ")) %>% 
          right_join(mz_selection)
        
        chrom_data_working <- chrom_data %>% filter(injection_id1 %in% inj_sel) %>% 
          group_by(injection,id1,injection_id1,scan_type,rt) %>% summarize(tic = sum(intensity)) %>% 
          ungroup() %>% 
          mutate(scan_type = as.character(factor(scan_type,levels=c("full","sim"),labels=c("Full Scan","SIM Scan"),ordered=T))) %>% 
          mutate(injection_label = paste0(injection,"\n",id1))
        
        chrom_data_downsampled <- chrom_data_working[seq(1,nrow(chrom_data_working),by=5),] %>% 
          group_by(scan_type) %>% 
          mutate(max_resp = max(tic))
        
        chromatogram_plot <- chrom_data_working %>%
          ggplot(aes(x=rt,y=tic,color=injection_label,customdata = scan_type)) +
          geom_point(data=chrom_data_downsampled,alpha=0) +
          geom_path() +
          scale_x_continuous(limits = c(min_rt,max_rt),expand = expand_scale(mult=0,add=0),n.breaks=25) +
          scale_y_continuous(expand = expand_scale(mult = c(0,1.1),add=0)) +
          facet_wrap(~scan_type,scales="free_y",ncol=1,strip.position = "right") +
          theme(legend.position = "top",
                strip.background = element_blank(),
                panel.border = element_rect(fill=NA,color="black")
          ) +
          labs(x="Retention Time (minutes)",
               y="Instrument Response (arbitrary units)",
               color="Injection")

        list("chromatogram_plot" = chromatogram_plot,
             "passcheck" = T)
      } else {
        list("chromatogram_plot" = ggplot()+annotate("text",x=0,y=0,label="No data loaded.",color="red")+theme_void(),
             "passcheck" = F)
      }
    }
    
    spectra_function <- function(rawdata,
                                 standard_info,
                                 gcms_names,
                                 inj_sel,
                                 rt_selection){
      
      if(!is.null(rt_selection)){
        spectral_selection <- rt_selection %>% 
          group_by(scan_type) %>% 
          summarize(lower = min(x),
                    upper = max(x)) %>% 
          ungroup() %>% 
          mutate(scan_type = ifelse(scan_type == "Full Scan","full","sim"))

        standard_spectral_data <- standard_info %>% select(line,injection) %>% distinct() %>%
          left_join(gcms_names,by = join_by(line, injection)) %>%
          left_join(rawdata,by = join_by(injection)) %>% 
          right_join(spectral_selection,by="scan_type") %>% 
          filter(rt >= lower & rt <= upper) %>%
          group_by(injection,mz) %>%
          summarize(total_intensity = sum(intensity)) %>% 
          group_by(injection) %>%
          mutate(scaled_intensity = total_intensity / max(total_intensity)) %>%
          group_by(mz) %>%
          summarize(mean_intensity = mean(scaled_intensity)) %>% 
          arrange(mean_intensity)
        
        sample_spectral_data <- data.frame(injection_id1 = inj_sel) %>%
          separate_wider_delim(cols = injection_id1,delim = " | ",
                               names=c("injection","id1")) %>%
          mutate(injection_id1 = paste0(injection,"_",id1)) %>% 
          left_join(rawdata,by="injection") %>%
          left_join(spectral_selection,by="scan_type") %>%
          filter(rt >= lower & rt <= upper) %>% 
          group_by(injection_id1,injection,mz) %>%
          summarize(total_intensity = sum(intensity)) %>%
          group_by(injection_id1,injection) %>%
          mutate(scaled_intensity = total_intensity / max(total_intensity)) %>% 
          arrange(-scaled_intensity) %>% 
          gather(intensity_type,value,c("total_intensity","scaled_intensity"))

        # print(standard_spectral_data,n=10)
        # print(sample_spectral_data,n=1000)
        
        sample_plot <- sample_spectral_data %>% 
          filter(intensity_type == "total_intensity") %>% 
          ggplot(aes(x=mz,y=value,color=injection_id1)) +
          geom_segment(aes(xend=mz,yend=0),position = position_dodge(width=1/10))
        
        list("sample_plot" = sample_plot,
             "passcheck" = T)
      } else {
        list("sample_plot" = ggplot()+annotate("text",x=0,y=0,label="No data loaded.",color="red")+theme_void(),
             "passcheck" = F)
      }
    }


    
    injection_list <- reactive({
      if(!is.null(ingest()$gcms_names)){
        ingest()$gcms_names %>%
          unite(injections,c(injection,id1),sep=" | ") %>% 
          .$injections} else return(NULL)
    })
    observe({
      updateCheckboxGroupInput(session,"injection_selection",choices = injection_list(),selected = injection_list()[1])
      updateNumericInput(session,"explorer_min_rt",value=ingest()$parameters$min_rt[1])
      updateNumericInput(session,"explorer_max_rt",value=ingest()$parameters$max_rt[1])
      updateTextInput(session,"explorer_fs_mz_range",value=ingest()$scan_info$mz_string[1])
      updateTextInput(session,"explorer_sim_mz_range",value=ingest()$scan_info$mz_string[2])
    })
    output$mz_syntax_description <- renderText({
      "Syntax:\nLeave empty to plot TIC. Enter ranges as numbers separated by a minus sign (e.g., 69 - 71),
      single values, or single values separated by commas (e.g., 71, 77, 78).
      You may also provide comma-separated sets of ranges (e.g., 71 - 77, 80, 90 - 400)."
    })
    
    chromatogram_explorer <- reactive({chromatogram_function(rawdata = ingest()$rawdata,
                                                             scan_info = ingest()$scan_info,
                                                             standard_info = ingest()$standard_info,
                                                             gcms_names = ingest()$gcms_names,
                                                             inj_sel = input$injection_selection,
                                                             min_rt = input$explorer_min_rt,
                                                             max_rt = input$explorer_max_rt,
                                                             full_mz = input$explorer_fs_mz_range,
                                                             sim_mz = input$explorer_sim_mz_range)})

    output$explorer_chromatogram_plot <- renderPlotly({
      plot <- chromatogram_explorer()$chromatogram_plot
      plotly_obj <- ggplotly(plot, dynamicTicks = T)
      
      gg <- plotly_obj %>% 
        event_register(event = "plotly_selected") %>%
        config(modeBarButtons = list(list("zoom2d", "pan2d", "select2d")),
               displaylogo = FALSE)
      if(chromatogram_explorer()$passcheck){
        # Extract correct facet-axis mapping from ggplotly()
        axis_map <- gg$x$data %>%
          purrr::map_dfr(~ data.frame(
            scan_type = .x$customdata[[1]],  # Extract correct scan_type from customdata
            xaxis = .x$xaxis %||% "x",
            yaxis = .x$yaxis %||% "y",
            max_y = max(.x$y)
          )) %>%
          distinct() %>% 
          group_by(scan_type,xaxis,yaxis) %>% 
          summarize(max_y = max(max_y))
        # Store mapping for use in selections
        session$userData$axis_map <- axis_map
        session$userData$plot_data <- gg$x$data
      }
      gg
    })
    
    output$explorer_spectra_range <- renderPrint({
      cat("No selection made. Use the 'Box Select' mode in the upper-right of plot to select range to display mass spectra. 
          You must draw the box over the plotted points on the chromatogram.")
    })

    observeEvent(event_data("plotly_selected"),{
      if(length(event_data("plotly_selected")) != 0){
        selected <- event_data("plotly_selected") %>% 
          rename(scan_type = customdata) %>%
          left_join(session$userData$axis_map,by="scan_type") %>%
          filter(x == min(x) | x == max(x))
        # print(selected)
        
        shape <- list(
          type = "rect",
          x0 = min(selected$x),
          x1 = max(selected$x),
          y0 = 0,
          y1 = selected$max_y[1],  # Explicitly taking the first value
          xref = selected$xaxis[1],
          yref = selected$yaxis[1],
          fillcolor = "rgba(255, 0, 0, 0.3)",  # Red with transparency
          line = list(width = 0)  # No border
        )
        
        plotlyProxy("explorer_chromatogram_plot", session) %>%
          plotlyProxyInvoke("relayout","shapes",list(shape))
        
        session$userData$explorer_selection <- selected
      } else {session$userData$explorer_selection <- NULL}
      
      output$explorer_spectra_range <- renderPrint({
        if (is.null(session$userData$explorer_selection)) {
          selected_output <- "No selection made. Use the 'Box Select' mode in the upper-right of plot to select range to display mass spectra.
        You must draw the box over the plotted points on the chromatogram."
        } else {
          x_range <- range(session$userData$explorer_selection$x)  # Extract min & max x-values
          selected_output <- paste("Selected RT Range:", round(x_range[1],2), "to", round(x_range[2],2), "minutes")
        }
        cat(selected_output)
      })
      
      session$userData$spectra_explorer <- spectra_function(rawdata = ingest()$rawdata,
                                                     standard_info = ingest()$standard_info,
                                                     gcms_names = ingest()$gcms_names,
                                                     inj_sel = input$injection_selection,
                                                     rt_selection = session$userData$explorer_selection)
      
      output$sample_spectra_plot <- renderPlotly({
        ggplotly(session$userData$spectra_explorer$sample_plot)
      })
      
    })
    

      
  
    

    


    
    
    

  
  }
  
}

shinyApp(ui, server = server)

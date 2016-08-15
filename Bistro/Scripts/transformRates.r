## R script to transform rates from bistro (in the log file)
## to MrBayes
## Usage:
## source("transformRates.r")
## > rates(logfile)
## Claudia August 2016

rates = function(file){
    lin <- readLines(file)
    i = 0
    for(i in 1:length(lin)){
        if(grepl("symmetric acceptance", lin[i]))
            break
    }
    p <- as.numeric(strsplit(lin[i+1],"\\s+")[[1]])
    r <- as.numeric(strsplit(lin[i+2],"\\s+")[[1]])
    p <- p[!is.na(p)]
    r <- r[!is.na(r)]
    den <- c(p[1]*p[2],p[1]*p[3],p[1]*p[4],p[2]*p[3], p[2]*p[4],p[3]*p[4])
    r <- r/den
    r <- r/r[6]
    print(paste0("prset statefreqpr=fixed(",p[1],",",p[2],",",p[3],",",p[4],");"))
    print(paste0("prset revmatpr=fixed(",r[1],",",r[2],",",r[3],",",r[4],",",r[5],",",r[6],");"))
}



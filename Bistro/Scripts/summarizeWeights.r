## R function to automatically summarize the weights from the table after bistro
## Usage: source("summarizeWeights.r")
## > summarizeWeights(file)
## Claudia August 2016

summarizeWeights = function(file){
    dat = read.table(file,header=TRUE)
    logweight = dat$logWt
    logw = logweight - max(logweight)
    w = exp(logw) / sum(exp(logw))
    dat$w <- w

    means=sort(with(dat, round(sapply( split( w,tree ),sum ),4) ), decreasing=TRUE)
    print("Posterior probabilities")
    if(length(means)>10){
        print(means[1:10])
        print(paste0("Only showing 10 trees out of ",length(means)))
    }else{
        print(means)
    }
    ses = rep(NA, length(means))
    for(i in 1:length(means)){
        se = with(dat, sqrt(sum(w^2*((tree==names(means)[i])-means[i])^2)))
        ses[i] <- se
    }
    ses = round(ses,4)
    print("Standard errors")
    sesnam = setNames(ses, names(means))
    if(length(sesnam)>10){
        print(sesnam[1:10])
        print(paste0("Only showing 10 trees out of ",length(sesnam)))
    }else{
        print(sesnam)
    }
}

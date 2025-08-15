# Get CPI data from FRED - there seems to be no function in Matlab that does an 
# FRED API call, so R it has to be *sigh*

rm(list=ls()) # clear memory

library(fredr)
library(tibble)

# category_id 9 are the CPI series, now get the childres
cats = fredr_category_children(category_id = 9)$id
ns = 0

for (i in 1:length(cats)){
  
  serNos = fredr_category_series(cats[i],
            filter_variable = "frequency",
            filter_value = "Monthly"
                                )

  for (j in 1:dim(serNos)[1]){

    # get only the seasonally adjusted series and those that are indices
    if (((serNos$seasonal_adjustment_short[j] == "SA") + (substr(serNos$units_short[j],1,5) == "Index")) == 2){
      ns = ns+1
      print(paste0("getting series ", j))
      series = fredr(serNos$id[j])
      write.csv(series,file=paste0("~/Dropbox (Personal)/Panel Forecasting/Data/CPI2023/series/", series[[1,2]],".csv"))
    }
  }
}
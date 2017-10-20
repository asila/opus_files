##.................................................................................
# Begin by sourcing where the read.opus file is stored

source('~/studies/Central_NRB_training/scripts/read.opus.R', chdir = TRUE)

#Change path to point to the folder with OPUS files

opus.files.directory <- "~/Dropbox/temp/share/soil_helth_status/OPUS_files"

#Set working directory to save the flat table to be created after conversions

setwd("~/studies/Rolf_Sommer")

#Give name to be used for saving the created table
  
file.name <- "Raw spectra.csv"

lst <- as.list(list.files(path=opus.files.directory, pattern=".[0-9]$", full.names=TRUE))

length(lst)

spectra<-c()

for ( i in 1:length(lst)){

spec <- opus(lst[[i]], speclib="ICRAF",plot.spectra=TRUE)

spectra<-rbind(spectra,spec)

}

#View part of spectra
spectra[1:5,1:8]
  
  
write.table(spectra,file=file.name,sep=",",row.names=FALSE)
#Read the file with eplots as a csv file

# Install required packages if not already done so.


is.installed<-function(anypkg){
	is.element(anypkg, installed.packages()[,1])
}

if (!is.installed(c("qrcode", "png","stringi"))){
	install.packages(c("qrcode", "png","stringi"))
}	

library(qrcode)
library(png)
library(grid)
library(stringi)

# Create folder where QR codes will be stored
dir.create("./QRcodes_images")
setwd("./QRcodes_images")

	#Read the file created from shiny app with alpha-numeric random codes.
	#change path and file name 	accordingly
	#qrc<-read.csv("d:/shiny app/new qr codes.csv") 
	qrc<-read.csv("~/new_qrcodes.csv") 

		qrcssn<-qrc[,1]
	qrcssn[1]
	n <-length(qrcssn)
	qrcodes<-c(1:n)
	for (i in 1:n){
    	x<-qrcssn[i]
    	y<-paste(x,sep="")
    	qrcodes[i]<-paste(y,collapse="",sep="")
    	fnam<-paste(qrcssn[i],".png",sep="")
    	png(file=fnam)
    	qrcode_gen(qrcodes[i],dataOutput=FALSE,plotQRcode=TRUE)
    	dev.off()
    	}

# Add the alpha numeric string to the image
	ll <- list.files(patt='*.png')
	imgs <- lapply(ll,function(x){
    img <- as.raster(readPNG(x))
    x.name <- gsub('(.*).png','\\1',x)
    png(file =paste(x.name,'.png',sep=''))
    grid.raster(img)
    grid.text(label = x.name,x=0.5,y=0.08,gp=gpar(cex=5))
    dev.off()
    })
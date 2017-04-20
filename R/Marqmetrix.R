## marqmetrix.R - import filter for Marqmetrix Raman specter files
## Author: Olaf Trygve Berglihn <oberg@sintef.no>
## Date: 2017-02-07

scan.txt.marqmetrix <- function(files = "*.txt", ..., label = list(),
				abscissa="Wavelength",
				ordinate="Counts") {

    ## set some defaults
    long <- list(files = files, ..., label = label)
    label <- modifyList(list(.wavelength = expression(Wavelength),
        spc = expression (Count)), label)

    ## find the files
    files <- Sys.glob(files)

    if (length(files) == 0) {
        warning("No files found.")
        return(new("hyperSpec"))
    }
    
    if (length(files) == 0) {
        warning("No files found.")
        return(new("hyperSpec"))
    }

    ## read the first file
    buffer  <- read.table(files[1], skip=14, header=TRUE, fill=TRUE)
    wavelength <- buffer[, abscissa]
    data <- scan.txt.marqmetrix.header(files[1])
    darkfile <- paste(dirname(files[1], data$processing, sep=""))
    darkbuffer <- read.table(darkfile, skip=14, header=TRUE, fill=TRUE)
    spc <- matrix(ncol = nrow(buffer), nrow = length(files))
    spc[1, ] <- buffer[, ordinate] - darkbuffer[, ordinate]

    ## read remaining files
    for (f in seq(along=files)[-1]) {
        buffer  <- read.table(files[f], skip=14, header=TRUE, fill=TRUE)
	hdr <- scan.txt.marqmetrix.header(files[f])
        darkfile <- paste(dirname(files[1], data$processing, sep=""))
        darkbuffer <- read.table(darkfile, skip=14, header=TRUE, fill=TRUE)

 	## Check wether they have the same wavelength axis
        if (! all.equal(buffer[, 1], wavelength))
            stop(paste(files[f], "has different wavelength axis."))
    
        spc[f, ] <- buffer[, ordinate] - darkbuffer[, ordinate]
        data <- rbind(data, hdr)
    }
    
    ## make the hyperSpec object
    new("hyperSpec", wavelength = wavelength, spc = spc, data = data,
        label = label)
}


scan.txt.marqmetrix.header <- function(file) {
    hdr <- scan(file, what ="raw", sep='\t', nlines=13)
    file <- file
    software.version <- hdr[min(which(match(hdr, "Software Version:")
				      == TRUE))+1]
    spectrometer.sn <- hdr[min(which(match(hdr, "Spectrometer SN:")
				     == TRUE))+1]
    date <- strptime(hdr[min(which(match(hdr, "Date:")==TRUE))+1],
    		  format="%m/%d/%Y %H:%M:%S %p", tz="CET")
    integration.time <- as.numeric(hdr[min(which(match(hdr, 
        "Integration time:")==TRUE))+1])
    averages <- as.numeric(hdr[min(which(match(hdr, "Averages:")==TRUE))+1])
    delay <- as.numeric(hdr[min(which(match(hdr, "Delay:") == TRUE))+1])
    sequence.number <- as.numeric(hdr[min(which(match(hdr, "Sequence Number:")
					       	== TRUE))+1])
    processing <- hdr[min(which(match(hdr, "Processing:") == TRUE))+1]
    background <- hdr[min(which(match(hdr, "Background:") == TRUE))+1]
    spectral.points <- as.numeric(hdr[min(which(match(hdr, "Spectral Points:")
					       	== TRUE))+1])
    laser.power <- hdr[min(which(match(hdr, "Laser Power:") == TRUE))+1]
    system.temp <- hdr[min(which(match(hdr, "System Temp:") == TRUE))+1]
    tec.temp <- hdr[min(which(match(hdr, "TEC Temp:") == TRUE))+1]	                    
    return(data.frame(file, software.version, spectrometer.sn, date,
        integration.time, averages, delay, sequence.number, processing,
        background, spectral.points, laser.power, system.temp, tec.temp))
}		       	 


## marqmetrix.R - import filter for Marqmetrix Raman specter files
## Author: Olaf Trygve Berglihn <oberg@sintef.no>
## Date: 2018-12-18

scan.txt.marqmetrix <- function(files = "*.txt", ..., label = list(),
				abscissa="RamanShift",
				ordinate="DarkSubtracted",
				darks=FALSE) {

    ## set some defaults
    long <- list(files = files, ..., label = label)
    label <- modifyList(list(.wavelength = expression(Raman ~ shift ~ "["*cm^{-1}*"]"),
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

    ## detect the number of lines for the header
    filebuf <- readLines(files[1])
    headlines <- which(! nzchar(filebuf))[1]

    ## Find the column titles, remove spaces
    colnames <- array(unlist(lapply(strsplit(filebuf[headlines+1], '\t'), gsub,
                       pattern = " ", replacement = "", fixed=TRUE)))
    
    ## read the first file
    buffer  <- read.table(files[1], skip=headlines+2, header=FALSE, fill=TRUE, col.names=colnames)
    wavelength <- buffer[, abscissa]
    data <- scan.txt.marqmetrix.header(files[1], headlines)

    spc <- matrix(ncol = nrow(buffer), nrow = length(files))
    
    if (darks && !identical(df$Processing, "dark subtract")) {
        darkfile <- strsplit(file.path(dirname(files[1]), data$Processing), ',')[[1]][1]
        darkbuffer <- read.table(darkfile, skip=headlines+1, header=TRUE, fill=TRUE)
	    spc[1, ] <- buffer[, ordinate] - darkbuffer[, ordinate]
    } else {
        spc[1, ] <- buffer[, ordinate]
    }

    ## read remaining files
    for (f in seq(along=files)[-1]) {
        buffer  <- read.table(files[f], skip=headlines+2, header=FALSE, fill=TRUE, col.names=colnames)
	hdr <- scan.txt.marqmetrix.header(files[f], headlines)

 	## Check wether they have the same wavelength axis
        if (! all.equal(buffer[, abscissa], wavelength))
            stop(paste(files[f], "has different wavelength axis."))

        if (darks && !identical(df$Processing, "dark subtract")) {
            darkfile <- strsplit(file.path(dirname(files[f]), hdr$Processing), ',')[[1]][1]
	    darkbuffer <- read.table(darkfile, skip=headlines+1, header=TRUE, fill=TRUE)
	    spc[f, ] <- buffer[, ordinate] - darkbuffer[, ordinate]
	} else {
	    spc[f, ] <- buffer[, ordinate]
	}

        data <- rbind(data, hdr)
    }
    
    ## make the hyperSpec object
    new("hyperSpec", wavelength = wavelength, spc = spc, data = data,
        label = label)
}


scan.txt.marqmetrix.header <- function(file, nlines) {
    hdr <- scan(file, what ="raw", sep='\t', nlines=nlines)

    keywords <- unlist(lapply(hdr[seq(1,length(hdr), by=2)],
                              function(x) {make.names(str_remove(x, ':$'))}))
    keywords <- c(keywords, "File")
    
    values <- hdr[seq(2,length(hdr),by=2)]
    values <- c(values, file)

    df <- data.frame(matrix(ncol=length(keywords), nrow=0))
    colnames(df) <- keywords
    df[1,] <- values
    return(df)
}		       	 


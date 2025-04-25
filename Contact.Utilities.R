# Utility Scripts to analyze hyvolution contact data

# Hyvolution matlab data output
    # 1. Stats csv
        # 2 variables: Consecutive time (contact event length), Image name (cell from which contact event from)
    # 2. Summary csv

data.to.read <- function(data.file) {
    #' data.to.read
    #'
    #' Receives an argument indicating which file to read, puts the file's data into a data frame and returns the data
    #'
    #' @param data.file a csv file containing contact raw matlab contact data (columns; ConsecutiveTime, ImageName. Each row is a contact event)
    #'
    #' @return data as a data frame

    # Defensive program
    # Stop and print "No csv file in directory" if no csv file exists
    if(file.exists(data.file) == FALSE) {
        stop("No csv file in directory")
    }

        # Convert external data to data frame
        contact.df <- read.csv(data.file)

        # Return data
        return(contact.df)
}

DwellTime.Freq.Pct.per.cell <- function(raw.contact.df, time) {
    #' DwellTimes.percent.per.cell
    #'
    #' Function calculates the absolute number and percentage of Dwell Times per cell
    #'
    #' @param contact.df data frame returned by contact.data.to.read
    #'
    #' @return pct_contact.csv containing columns; ConsecutiveTime, ImageName, Freq (absolute contact frequency), Pct (% contact frequency). Each row is for 1 DwellTime for a particular cell.

    # Create framework dataframe
    # Important so that all timpoints have value even if 0 for all images
    framework.df <- data.frame(ConsecutiveTime = 1:time, ImageName = "framework")

    # Add framework dataframe to raw.contact.df
    merged.df <- rbind(raw.contact.df, framework.df)

    # Create ftable of raw contact duration frequencies from merged.df
    ftable.table <- ftable(merged.df)

    # Convert ftable to data frame
    ftable.df <- as.data.frame(ftable.table)

    # Remove framework rows
    ftable.df <- subset(ftable.df, ImageName!="framework")
    ftable.df$ImageName <- factor(ftable.df$ImageName)

    # Make ConsecutiveTime column integer
    ConsecutiveTime <- as.integer(ftable.df$ConsecutiveTime)
    ftable.df$ConsecutiveTime <- ConsecutiveTime

    # Slice data frame to ImageName column
    ImageNames <- ftable.df$ImageName

    # Identify the different Images
    Images <- unique(ImageNames)

    #Initialize cumulative list
    composite.list <- list()

    # Apply subsequent loop to each image
    for(each.Image in Images) {

        # Slice ftable.df to rows of Image being looped
        each.Image.df <- ftable.df[ftable.df$ImageName == each.Image ,]

        # Define function to calculate percent frequency for Image being looped
        pctfn <- function(x) {x/sum(each.Image.df$Freq)*100}

        # Apply pctfn to Freq column corresponding to Image being looped
        Pct <- tapply(each.Image.df$Freq, each.Image.df$ConsecutiveTime, pctfn)

        # Convert tapply list output to a dataframe
        pct.df <- as.data.frame(Pct)

        # Add data to list
        composite.list[[each.Image]] <- pct.df
    }

    # Combine listed outputs into one dataframe by adding rows together
    # Yeilds data frame containing 1 column of % contact with row name ImageName.ConsecutibeTime
    composite.df <- do.call(rbind, composite.list)


    # Combine composite.df (containing % contact calculations) with ftable.df (containing ImageName and ConsecutiveTime columns)
    # Because row order is same % contact will match with appropriate ImageName and ConsecutiveTime columns
    absolute.pct.df <- cbind(ftable.df, composite.df)

    write.csv(absolute.pct.df, "absolute_pct_contact.csv")

    return(absolute.pct.df)
}

Freq.per.cell <- function(absolute.pct.df) {
    #' DwellTimes.percent.per.cell
    #'
    #' Function creates output for curve fitting. Freqency distributions for each cell are put into own rows with ImageName as Column name
    #'
    #' @param absolute.pct.df data frame returned by DwellTime.Freq.Pct.per.cell
    #'
    #' @return Freq_per_cell.csv containing Freq columns per cell

    # Split Freq column based on ImageName. Output is a list
    Freq.per.cell.list <- split(absolute.pct.df$Freq, absolute.pct.df$ImageName)

    # Combine list as matrix
    Freq.per.cell.matrix <- sapply(Freq.per.cell.list, cbind)

    # Convert matrix to data frame
    Freq.per.cell.df <- as.data.frame(Freq.per.cell.matrix)

    # Export data frame as a csv
    write.csv(Freq.per.cell.df, "Freq_per_cell.csv")
}

XBinned.PctPO.per.cell <- function(summary.df, absolute.pct.df, x) {
    #' FINAL.PctPO.per.cell
    #'
    #' Function calculates the percentage of peroxisomes in contact for = or > x timepoints
    #'
    #' @param summary.df data frame returned by contact.data.to.read
    #' @param absolute.pct.df data frame returned by DwellTime.Freq.Pct.per.cell
    #'
    #' @return pct_contact.csv containing columns; ConsecutiveTime, ImageName, Freq (absolute contact frequency), Pct (% contact frequency)


    # Slice summary.df to tracked_samples column (gives number of peroxisomes imaged)
    total.num.peroxisomes <- summary.df$tracked_samples

    # Slice absolute.pct.df to rows where ConsectutiveTime is >= 4
    x.freq.pct.df <- absolute.pct.df[absolute.pct.df$ConsecutiveTime >= x, ]

    # Calculate sum frequency when ConsecutiveTime >= 4 per image and convert result to data frame
    sum.Freq.xUp <- tapply(x.freq.pct.df$Freq, x.freq.pct.df$ImageName, sum)
    sum.Freq.xUp.df <- as.data.frame(sum.Freq.xUp)


    # Normalize to number of PO
    # Because image order is same peroxisome.num will match with appropriate ImageName
    pctPO <- (sum.Freq.xUp.df$sum.Freq.xUp/total.num.peroxisomes)*100

    Final.xUp.df <- cbind(sum.Freq.xUp.df, total.num.peroxisomes, pctPO)

    write.csv(Final.xUp.df, paste0("%PO_", x, "up.csv"))

    return(Final.xUp.df)
}

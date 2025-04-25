# Input data: Hyvolution matlab data CSV output
    # 1. Stats csv
        # 2 variables: Consecutive time (contact event length), Image name (cell from which contact event from)
    # 2. Summary csv


# Defensive program
# Stop if no csv file exists
if(length(dir(pattern = ".csv")) == 0) {
      stop("No csv file in directory")
}


# File containing the functions
source("Contact.Utilities.R")


# Create a data frame for the raw data and return the data assigned to raw.contact.df
for (filename in dir(pattern = "stats")) {
    raw.contact.df <- data.to.read(filename)
}

# Create a data frame for the raw data and return the data assigned to summary.df
for (filename in dir(pattern = "summary")) {
    summary.df <- data.to.read(filename)
}


# Calculate the absolute number and percentage of contact events per dwell time (where 10 is last timepoint) per cell and export as a csv
absolute.pct.df <- DwellTime.Freq.Pct.per.cell(raw.contact.df, 10)

# Create output for curve fitting. Freqency distributions for each cell are put into own rows with ImageName as column name
Freq.per.cell(absolute.pct.df)

# Calculate the percentage of peroxisomes in contact for >=6 timepoints per cell and export as a csv
XBinned.PctPO.per.cell(summary.df, absolute.pct.df, 6)

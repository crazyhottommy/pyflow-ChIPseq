#!/usr/bin/env Rscript

# R script for downloading files from SRA using ascp.
# See help for usage. 
# this script is from https://github.com/lazappi/binf-scripts/blob/master/sraDownload.R thanks Luke Zappia for sharing!

library("optparse")

# Setup options
opt.list = list(
    make_option(c("-d", "--database"),
                type    = "character",
                default = "./SRAmetadb.sqlite", 
                help    = paste("Path to SRA database file, if does not exist",
                                "database will be download to this location.",
                                "Default is 'SRAmetadb.sqlite'.")),
    make_option(c("-o", "--outpath"),
                type    = "character",
                default = getwd(),
                help    = paste("Path to save downloaded SRA files.",
                                "Default is CWD.")),
    make_option(c("-a", "--ascpCMD"),
                type    = "character",
                default = "ascp -QT -l 300m -i /usr/local/aspera/connect/etc/asperaweb_id_dsa.openssh",
                help    = paste("Command for running ascp.",
                                "Default is 'ascp -QT -l 300m -i",
                                "/usr/local/aspera/connect/etc/asperaweb_id_dsa.openssh'.")),
    make_option(c("-t", "--type"),
                type    = "character",
                default = "sra",
                help    = paste("Type of files to download. Options are",
                                "'sra' (default) or `fastq`."))
)

# Read options from command line
parser     <- OptionParser(usage = "%prog [options] SRA1, SRA2, ...",
                            option_list = opt.list)
arguments  <- parse_args(parser, positional_arguments = c(1, Inf))
opt        <- arguments$options
ids        <- arguments$args

# Load other libraries here so it only happens if options are good
suppressPackageStartupMessages(library("SRAdb"))
library("DBI")

# Check if SRA databases exists and download if needed

message("Checking if database exists...")
if (!file.exists(opt$database)) {
    message("Database file not found.")
    db.path  <- strsplit(opt$database, split = "/")[[1]]
    db.file  <- db.path[length(db.path)]
    db.path  <- do.call(file.path, as.list(db.path[-length(db.path)]))
    if (!dir.exists(db.path)) {
        dir.create(db.path)
    }
    message(paste("Downloading database to", file.path(db.path, db.file)))
    sql.file <- getSRAdbFile(destdir  = db.path,
                             destfile = paste0(db.file, ".gz"))
} else {
    message(paste("Database file found at", opt$database))
    sql.file <- opt$database
}

# Connect to the database
message("Connecting to database...")
sra.con <- dbConnect(RSQLite::SQLite(), sql.file)

# Create output directory if it doesn't exist
message("Checking if output directory exists...")
if (!dir.exists(opt$outpath)) {
    message(paste("Creating output directory at", opt$outpath))
    dir.create(opt$outpath)
}

# For each SRA id download the metadata and zipped fastqs
message("Processing SRA IDs...")
num.ids <- length(ids)
for (idx in 1:num.ids) {
    id <- ids[idx]
    message(paste0("Processing ", id, ", ", idx, " of ", num.ids, "..."))
    meta <- listSRAfile(id, sra.con, fileType = opt$type)
    study <- unique(meta$study)
    query <- paste("SELECT * FROM experiment WHERE study_accession IN (\'",
                   study,"\')",sep="")
    message(paste("Query:", query))
    exp <- dbGetQuery(sra.con, query)
    meta <- merge(meta, exp, by.x="experiment", by.y="experiment_accession")
    # Clean up sample titles so they can be used as file names
    meta$title <- gsub(" ","_",gsub(":","",meta$title))
    message(paste(nrow(meta), "SRA records found"))
    meta.path <- file.path(opt$outpath, paste0(id, ".metadata"))
    write.table(meta,
                file      = meta.path,
                quote     = FALSE,
                sep       = "\t",
                row.names = FALSE)
    message(paste("Metadata saved to", meta.path))
    message("Downloading files...")
    getSRAfile(id,
                sra.con,
                destDir  = opt$outpath,
                fileType = opt$type,
                srcType  = "fasp",
                ascpCMD  = opt$ascpCMD)
    message("Renaming files...")

    files <- basename(meta$ftp)

#    file.rename(file.path(opt$outpath, files),
#               file.path(opt$outpath, paste(meta$title, files, sep = ".")))
    message(paste(id, "finished"))
}

message("Done!")
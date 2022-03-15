makeblastdb <- function(infile, 
                        input_type="fasta", 
                        dbtype="nucl", 
                        hash_index = FALSE,
                        workdir = "."){
  
  if (!file.exists(infile)) stop("infile does not exists.")
  infile <- normalizePath(infile)
  
  workdir <- normalizePath(workdir)
  wd <- getwd()
  setwd(workdir)
  on.exit(setwd(wd))
  
  out <- sub("[.]\\w+$", "", basename(infile))
  cmd <- paste(
    "makeblastdb",
    "-in", infile,
    "-out", out,
    "-input_type", input_type,
    "-dbtype", dbtype,
    ifelse(hash_index, "-hash_index", "")
  )
  system(cmd)
  ext <- c("ndb", "nhr", "nin", "not", "nsq", "ntf", "nto")
  if (hash_index){
    hix <- c("nhd", "nhi", "nog")
    ext <- c(ext, hix)
  }
  outfiles <- normalizePath(paste(out, ".", ext, sep = ""))
  fex <- file.exists(outfiles)
  # print(fex)
  if (!all(fex)){
    # print(fex)
    stop("Error creating blast database.")
  }
  return(c(file.path(workdir, out), outfiles))
}

blastn <- function(query, 
                   subject, 
                   task = c('blastn', 'blastn-short', 'dc-megablast', 'megablast' ,'rmblastn'),
                   evalue=1e-3, 
                   outfmt6 = "qseqid sseqid pident evalue gaps length qcovs",
                   num_threads = 1L,
                   workdir = ".",
                   ...){
  
  task <- match.arg(task, 
                    choices = c('blastn', 'blastn-short', 'dc-megablast', 'megablast' ,'rmblastn'), 
                    several.ok = FALSE)
  
  if (!file.exists(query)) stop("query does not exists.")
  if (!file.exists(subject)) stop("subject does not exists.")
  
  query <- normalizePath(query)
  subject <- normalizePath(subject)
  
  wd <- getwd()
  setwd(workdir)
  on.exit(setwd(wd))
  
  db <- try(makeblastdb(infile = subject, ...))
  if (class(db) != "try-error") on.exit(file.remove(db[-1]), add = TRUE)
  
  outfile <- paste0(
    "blastn_",
    sub("[.]\\w+$", "", basename(query)),
    "_vs_",
    basename(db[1]),
    ".tsv"
  )
  
  cmd <- paste(
    "blastn",
    "-task", task,
    "-query", query,
    "-db", db[1],
    "-evalue", evalue,
    "-outfmt", paste("'6 ", outfmt6, "'", sep = ""),
    "-out", outfile,
    "-num_threads", num_threads
  )
  
  message(cmd)
  system(cmd)
  
  if (file.info(outfile)$size != 0){
    cols <- strsplit(outfmt6, " ")[[1]]
    df <- read.csv(outfile, 
                   header = FALSE, 
                   sep = "\t", 
                   col.names = cols)
  }else{
    df <- data.frame()
  }
  attr(df, "blastn_outfile") <- outfile
  return(df)
}

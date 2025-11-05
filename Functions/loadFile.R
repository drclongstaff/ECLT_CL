#Find extension of file to load
getExtension <- function(file) {
  filename <- basename(file)
  last_dot <- regexpr("\\.([^.]+)$", filename)
  if (last_dot == -1) {
    return("")
  } else {
    return(substr(filename, last_dot + 1, nchar(filename)))
  }
} 
#Load text or csv file
load_file <- function(NAME, PATH, SHEET){
  ext <- getExtension(NAME)
  switch(ext,
         txt = read.delim(PATH),
         csv = read.csv(PATH),
         validate("Invalid file. Please upload a data file")
  )
}

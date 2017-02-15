

library(xlsx)
library(yaml)
library(XML)
library(RCurl)
library(rjson)
library(rcrossref)


rs  <- read.table("www/datasets.txt", sep="\t", header = T, stringsAsFactors = F)



# GET THE CITATION INFO FROM CROSSREF
for(i in 1:nrow(rs)){
  if(rs$doi[i] != "0" & rs$doi[i] != "_"){
    # Get more info with the pmid id
    tryCatch({
      # Get the XML data
      art <-cr_cn(dois = rs$doi[i], format = "citeproc-json")

      # # If journal does not exist, create the entry
      rs$journal[i] <- art$`container-title`
      
      # Fill the author table
      author_string <- ""
      fauthor = ""
      authors <- art$author
      for(k in 1:(length(authors$family))){
        lname <- authors$family[k]
        fname <- authors$given[k]
        
        fname <- strsplit(fname, " ")[[1]]
        forename <- NULL
        sep=""
        for(j in 1:length(fname)){
          if(j > 1) sep = " "
          forename <- paste(forename, substring(fname[j], 1, 1), sep=sep)
        }
        if(k == 1){
          author_string <- paste0(lname, " ", forename)
          rs$first_author[i] <- paste0(lname, " ", forename)
        } else author_string <- paste0(author_string, ", ", lname, " ", forename)
      }
      
      author_string <- gsub("<U+00E1>", "a", author_string)
      author_string <- gsub("<U+00C1>", "A", author_string)
      author_string <- gsub("<U+00E8>", "e", author_string)
      author_string <- gsub("<U+00E8>", "e", author_string)
      rs$author[i] <- author_string
      
      if(!is.null(art$title)) rs$title[i] <- art$title
      if(!is.null(art$volume)) rs$volume[i] <- art$volume
      if(!is.null(art$page)) rs$pages[i] <- art$page
      if(!is.null(art$`published-print`$`date-parts`[1])) rs$year[i] <- art$`published-print`$`date-parts`[1]
      
      message(paste("-----------  ",i,"  --------", rs$name[i]))
      
      # Store the failed PMIDS
    }, warning = function(w) {
      print(w)
      message(paste(rs$doi[i], " crossref failed: ", rs$name[i]))
      return(w)
    }, error = function(e) {
      print(e)
      message(paste(rs$doi[i], " crossref failed: ", rs$name[i]))
      return(e)
    })
  }
}
remove(art, author_string, authors, fauthor, fname, forename, i, j, k, lname, sep)

write.table(rs, "www/datasets.txt", quote=F, row.names = F, sep="\t")


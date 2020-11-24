## ---------------------------
##
## Script name: COVID-19 PubMed Analysis
##
## Purpose of script: Code to query PubMed for all articles related to COVID-19
##
## Authors: Saumya Gurbani, Mohleen Kang, Jordan Kempker
##
## Date: 11/24/2020
##
## Email: saumya.gurbani@emory.edu
##
## ---------------------------

# Libraries
library("rjson")
library("xml2")

# Clear workspace
remove(list = ls())


# Helper function to convert Entrez XML Date node to R Date format; the node has 3 subfields for Year, Month, and Day
EntrezXMLDate <- function(mydate) {
  return(as.Date(paste(mydate$Year,mydate$Month,mydate$Day,sep="-"),format="%Y-%m-%d"))
}

# Helper function to convert JournalIssue XML Date node to R Date format; the node has 3 subfields for Year, Month, and Day
JournalIssueXMLDate <- function(mydate) {
  return(as.Date(paste(mydate$Year,mydate$Month,mydate$Day,sep="-"),format="%Y-%b-%d"))
}

# Helper function to convert MeshHeading to a character list of Mesh terms
MeshList <- function(mesh_xml) {
  mesh = list()
  for(i in 1:length(mesh_xml)) {
    mesh <- rbind(mesh, as.character(unlist(mesh_xml[k]$MeshHeading$DescriptorName)))
  }
  
  return( mesh )
}

# Get list ----------------------------------------------------------------

# timing
ptm <- proc.time()

assign("last.warning", NULL, envir = baseenv())

# Variables
base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/" # Entrez API base
search_params = "esearch.fcgi?db=pubmed&retmode=json&retmax=100000&term="
fetch_params = "efetch.fcgi?db=pubmed&rettype=abstract&id="
batch_size = 200

#search query is the default COVID-19 query at https://www.nlm.nih.gov/index.html#Novel_Coronavirus
search_query = '((wuhan[All%20Fields]%20AND%20("coronavirus"[MeSH%20Terms]%20OR%20"coronavirus"[All%20Fields]))%20AND%202019/12[PDAT]%20:%202030[PDAT])%20OR%202019-nCoV[All%20Fields]%20OR%202019nCoV[All%20Fields]%20OR%20COVID-19[All%20Fields]%20OR%20SARS-CoV-2[All%20Fields]'


### Run query using Entrez esearch
print("Querying Entrez using search terms...")
query_url = paste(base_url, search_params, search_query, sep="")
search_results =  rjson::fromJSON(file=query_url)

num_results = as.integer(search_results$esearchresult$count)
print(paste("...done.", num_results, "results returned."))

# Get list of ID's
idlist = search_results$esearchresult$idlist


# Fetch -------------------------------------------------------------------


# prepare results frame
#data = data.frame(col.names=c("AD","AU","CRDT","DEP","DP","EDAT","GR","LA","MH","PHST","PL","PMC","PMID","PST","PT","TA","TYPE"))
data = list()

# now we need to query Entrez for each of these ID's using efetch
# this can occur in batches of ~200 at a time

ind = 0
while(ind < length(idlist)) {
  # get the next batch_size of IDs
  id_batch= paste(idlist[(ind+1):(ind+batch_size)],collapse=",")
  query = paste(base_url, fetch_params, id_batch, sep="")
  
  results = xml2::as_list( xml2::read_xml(query) )[1]$PubmedArticleSet
  
  print(paste("Starting index ", ind))
  
  for(j in 1:batch_size) {
    
    if( j > length(results) ) {
      break()
    }
    
    res = results[j]
    
    # create blank item
    item = list()
    
    # determine if this an Article or BookArticle
    if(attributes(res) == "PubmedBookArticle") {
      next()
    } else {
      res = res$PubmedArticle
      
      #item$TYPE = "article"
      
      # set the PubMed ID (PMID)
      item$PMID = as.character(res$MedlineCitation$PMID)
      
      # set Language (LA)
      item$LA = as.character(res$MedlineCitation$Article$Language)
      
      # set Article date of publication (DP)
      #item$DP = EntrezXMLDate(res$MedlineCitation$Article$ArticleDate)
      
      # set the Publication status (PST)
      item$PST = as.character(res$PubmedData$PublicationStatus)
      
      # set the list of publication types (PT), with a max of 2
      item$PT = unname(unlist(res$MedlineCitation$Article$PublicationTypeList))
      if(length(item$PT) > 1) {
        item$PT2 = item$PT[2]
        item$PT = item$PT[1]
      } else {
        item$PT2 = ""
      }
      
      # set the journal title abbreviation (TA)
      item$TA = as.character(res$MedlineCitation$MedlineJournalInfo$MedlineTA)
      
      # set the place (country) of publication (PL)
      item$PL = as.character(res$MedlineCitation$MedlineJournalInfo$Country)
      
      # set the Other Terms from the KeywordList (OT)
      item$OT = as.list(unname(unlist(res$MedlineCitation$KeywordList)))
      
      # set the GrantList, if it exists
      if(!is.na(match("GrantList", attributes(res$MedlineCitation$Article)$names))) {
        grants = list()
        for(k in 1:length(res$MedlineCitation$Article$GrantList)) {
          grant = list()
          grant_xml = res$MedlineCitation$Article$GrantList[k]$Grant
          grant$ID = as.character(unlist(grant_xml$GrantID))
          grant$Agency = as.character(unlist(grant_xml$Agency))
          
          grant$Institute = as.character(unlist(grant_xml$Institute))
          if(length(grant$Institute) == 0) {
            grant$Institute = NA
          }
          
          grant$Country = as.character(unlist(grant_xml$Country))
          if(length(grant$Country) == 0) {
            grant$Country = NA
          }
          
          grants <- rbind(grants, grant)
        }
        
        item$GR = grants
      } else {
        item$GR = NA
      }
    }
    
    # set the list of Mesh terms (MH), if exists
    if( !is.na(match("MeshHeadingList", attributes(res$MedlineCitation)$names)) ) {
      item$MH = MeshList(res$MedlineCitation$MeshHeadingList)
    } else {
      item$MH = NA
    }
    
    # Set the PMC identifier (PMC), if exists
    item$PMC = NA
    if( !is.na(match("ArticleIdList", attributes(res$PubmedData)$names)) ) {
      ail = res$PubmedData$ArticleIdList
      for(k in 1:length(ail)) {
        if(attr(ail[k]$ArticleId,"IdType") == "pmc") {
          item$PMC = as.character(unlist(ail[k]))
          break()
        }
      }
    }
    
    # Set the Entrez Date (EDAT), if exists
    item$EDAT = NA
    if( !is.na(match("History", attributes(res$PubmedData)$names)) ) {
      phist = res$PubmedData$History
      for(k in 1:length(phist)) {
        if(attr(phist[k]$PubMedPubDate,"PubStatus") == "entrez") {
          item$EDAT = as.character(EntrezXMLDate(phist[k]$PubMedPubDate))
          break()
        }
      }
    }
    
    # Set the Date of Electronic Publication (DEP), if exists
    item$DEP = NA
    if( !is.na(match("ArticleDate", attributes(res$MedlineCitation$Article)$names)) ) {
      item$DEP = as.character(EntrezXMLDate(res$MedlineCitation$Article$ArticleDate))
    }
    
    # Set the Journal Date of Publication (DP), if exists
    item$DP = NA
    if( !is.na(match("Journal", attributes(res$MedlineCitation$Article)$names)) ) {
      item$DP = as.character(JournalIssueXMLDate(res$MedlineCitation$Article$Journal$JournalIssue$PubDate))
    }
    
    # add to data frame
    data <- rbind(data, item)
  }
  
  ind = ind + batch_size
}

# Create a Data Frame from data
rownames(data) <- c()
pub_db <- as.data.frame(data)

# Create factors and remove listing within columns
pub_db$PMID <- unlist(pub_db$PMID)
pub_db$LA <- as.factor(unlist(pub_db$LA))
pub_db$PST <- as.factor(unlist(pub_db$PST))
pub_db$PT <- as.factor(as.character(pub_db$PT))
pub_db$PT2 <- as.factor(as.character(pub_db$PT2))

# Remove 0 length elements in PL and TA
pub_db$PL[lengths(pub_db$PL) == 0] <- NA
pub_db$PL <- as.factor(unlist(pub_db$PL))

pub_db$TA[lengths(pub_db$TA) == 0] <- NA
pub_db$TA <- as.factor(unlist(pub_db$TA))

# Create sortable Date columns
pub_db$DP[lengths(pub_db$DP) == 0] <- NA
pub_db$DP <- as.Date(unlist(pub_db$DP))

pub_db$DEP[lengths(pub_db$DEP) == 0] <- NA
pub_db$DEP <- as.Date(unlist(pub_db$DEP))

pub_db$EDAT[lengths(pub_db$EDAT) == 0] <- NA
pub_db$EDAT <- as.Date(unlist(pub_db$EDAT))

print(paste("Total time: ", round((proc.time() - ptm)[3],1), "(per record: ", round((proc.time()-ptm)[3]/ind,4), ")"))

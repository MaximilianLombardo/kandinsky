
#Open Targets - target query

endpoint <- "https://platform-api.opentargets.io/v3/platform/public/search?q=SNCA"

res <- jsonlite::fromJSON(endpoint)
res$data

#Open targets Private target endpoint

private.target.enpoint <- "https://platform-api.opentargets.io/v3/platform/private/target/ENSG00000145335"
res.private <- jsonlite::fromJSON(private.target.enpoint)


#Chembl api query

chembl.endpoint <- "https://ebi.ac.uk/chembl/api/data/target/"
    
res <- jsonlite::fromJSON(chembl.endpoint)

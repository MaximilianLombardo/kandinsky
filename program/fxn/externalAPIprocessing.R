#SCNA api query examples

#Open Targets - target query

endpoint <- "https://platform-api.opentargets.io/v3/platform/public/search?q=SNCA"

res <- jsonlite::fromJSON(endpoint)
res$data

#Open targets Private target endpoint

private.target.enpoint <- "https://platform-api.opentargets.io/v3/platform/private/target/ENSG00000145335"
res.private <- jsonlite::fromJSON(private.target.enpoint)

#Open targets - get chemical matter associated with a specific target
target.chembl.endpoint <- "https://api.opentargets.io/v3/platform/public/evidence/filter?target=ENSG00000145335&datasource=chembl"
res.chembl <- jsonlite::fromJSON(target.chembl.endpoint)

#Chembl api query



#Can't get this working?
chembl.endpoint <- "https://www.ebi.ac.uk/chembl/api/data/target/search?q=ENSG00000145335"
chembl.relation.endpoint <- "https://www.ebi.ac.uk/chembl/api/data/target_relation/CHEMBL6152"#CHEMBL6152
chembl.target <- jsonlite::fromJSON(chembl.endpoint)
chembl.relation <- jsonlite::fromJSON(chembl.relation.endpoint)

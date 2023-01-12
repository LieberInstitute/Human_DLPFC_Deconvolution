#!/usr/bin/env R

# Author: Sean Maden
#
# Inspect properties of recount3 datasets for DLPFC RO1 deconvolution project.
#
# 

library(recount3)
dir.name <- "datasets"

#----------
# load data
#----------
# get dataset
# note: csv obtained using "human", "dorsolateral" search terms
csv.fname <- "recount3_selection_2022-09-06 18_10_31.csv"
csv.fpath <- file.path(dir.name, csv.fname)
csv <- read.csv(csv.fpath)

#-----------------
# helper functions
#-----------------
map_libprep <- function(cdi, 
                        libprep.cnv = c("sra.library_construction_protocol", 
                                        "sra.study_abstract"), 
                        grepl.vectors = c(grepl.pattv, 
                                          grepl.pattv.polya, 
                                          grepl.pattv.ribozero), 
                        cnv = c("libprep.available", "pa.true", "rz.true")){
  # map_libprep
  # 
  # Gets the libprep condition outcome from srp coldata.
  # 
  # arguments:
  # cdi : coldata from rse
  # libprep.cnv : vector of variables, colnames in cdi to inspect
  # grepl.vectors : grepl pattern vectors corresponding to cnv items
  # cnv : new variable names corresponding to grepl.vectors items
  # 
  # returns:
  # dfi, data.frame of condition outcomes
  # 
  dfi <- do.call(cbind, lapply(grepl.vectors, function(gi){
    apply(do.call(cbind, lapply(libprep.cnv, function(ci){
      grepl(gi, cdi[,ci])
    })), 1, function(ri){
      ifelse(length(which(ri))>0, TRUE, FALSE)
    })
  }))
  colnames(dfi) <- cnv
  return(dfi)
}

#-------------
# get rse data
#-------------
lsubsets <- lapply(srpv, function(srpi){subset(human_projects, 
                                               project == srpi & 
                                                 project_type == "data_sources")})
lrse <- lapply(lsubsets, function(subi){create_rse(subi)})
names(lrse) <- srpv

# summary
for(ri in lrse){print(dim(ri))}
# [1] 63856    82
# [1] 63856    62
# [1] 63856    52
# [1] 63856    52
# [1] 63856    40
# [1] 63856    39
# [1] 63856    36
# [1] 63856    24
# [1] 63856    20
# [1] 63856    12
# [1] 63856     5

#-------------
# get metadata
#-------------
# coldata
lcd <- lapply(lrse, function(rsei){
  dfi <- as.data.frame(colData(rsei));dfi$srp <- srpi; dfi
})
names(lcd) <- srpv
# summary
for(ci in lcd){print(dim(ci))}
# [1]  82 176
# [1]  62 176
# [1]  52 176
# [1]  52 176
# [1]  40 176
# [1]  39 176
# [1]  36 169
# [1]  24 176
# [1]  20 169
# [1]  12 176
# [1]   5 176

# expanded sra attributes
lsr <- lapply(lrse, function(rsei){
  dfi <- as.data.frame(colData(expand_sra_attributes(rsei)))
  dfi$srp <- srpi; dfi
})
names(lsr) <- srpv
# summary
for(ci in lsr){print(dim(ci))}
# [1]  82 180
# [1]  62 179
# [1]  52 188
# [1]  52 181
# [1]  40 190
# [1]  39 181
# [1]  36 178
# [1]  24 180
# [1]  20 173
# [1]  12 182
# [1]   5 180
# view unique sra attribute variables
unique(as.character(unlist(lapply(lsr, function(ci){
  colnames(ci)[grepl("sra_attribute", colnames(ci))]
}))))
# [1] "sra_attribute.brain_region"         "sra_attribute.disease_stage"        "sra_attribute.source_name"         
# [4] "sra_attribute.tissue"               "sra_attribute.disease_state"        "sra_attribute.age"                 
# [7] "sra_attribute.brainid"              "sra_attribute.disease_status"       "sra_attribute.percentexonicmapping"
# [10] "sra_attribute.ph"                   "sra_attribute.race"                 "sra_attribute.rin"                 
# [13] "sra_attribute.Sex"                  "sra_attribute.totalnummappedreads"  "sra_attribute.totalnumreads"       
# [16] "sra_attribute.identifier"           "sra_attribute.biomaterial_provider" "sra_attribute.BioSampleModel"      
# [19] "sra_attribute.brain_number"         "sra_attribute.cell_subtype"         "sra_attribute.degradation_time"    
# [22] "sra_attribute.disease"              "sra_attribute.ethnicity"            "sra_attribute.flowcell"            
# [25] "sra_attribute.isolate"              "sra_attribute.library_type"         "sra_attribute.RIN"                 
# [28] "sra_attribute.sex"                  "sra_attribute.individual_id"        "sra_attribute.processing_date"     
# [31] "sra_attribute.case.control"         "sra_attribute.Cause_of_death"       "sra_attribute.manner_of_death"     
# [34] "sra_attribute.germinal_zone"        "sra_attribute.cell_type"            "sra_attribute.cell_line"           
# [37] "sra_attribute.differentiation_day"  "sra_attribute.gender"               "sra_attribute.genotype"            
# [40] "sra_attribute.pathology"

# save
# save coldata
cd.list.fname <- "list-coldata_dorsolateral-human_recount3.rda"
save(lcd, file = file.path(dir.name, cd.list.fname))
# save coldata csv
# cd.fname <- "coldata_dorsolateral-human_recount3.csv"
# save(do.call(rbind, lcd), file = cd.fname)

#---------------------------
# map sra_attribute metadata
#---------------------------
# set variable mappings
lvar <- list(cell_info = c("cell"),
             disease_info = c("disease","pathology","case"),
             sex = c("sex","gender"), tissue = c("tissue","region"),
             demographics = c("race","ethnicity"), age = c("age"),
             identifier = c("id", "ID"), death_details = c("death"),
             library = c("library"))

# get regex patterns
start.str <- ".*"; end.str <- ".*"
lgrep <- lapply(lvar, function(li){
  paste0(unlist(lapply(li, function(ii){
    charv <- unlist(strsplit(ii, "")); letter <- charv[1] # get first letter
    which.letter <- paste0("(", letter, "|", toupper(letter), ")") # 1st chr cap
    notfirst.str <- paste0(charv[2:length(charv)], collapse = "") # chr remain
    char.str <- paste0(which.letter, notfirst.str) # new str
    start.str <- ifelse(ii == "age", paste0(start.str, "(^|\\_|!s)"), start.str)
    paste0(start.str, char.str, end.str, collapse = "")
  })), collapse = "|")
})

# get mappings
# get colname mappings by srp
lmap <- lapply(lsr, function(sri){
  df.sri <- sri[,grepl("^sra_attribute.*",colnames(sri))]
  cnv <- colnames(df.sri)
  lmapi <- lapply(seq(length(lgrep)), function(ii){
    message(ii)
    varname <- names(lgrep)[ii]; patti <- lgrep[[ii]]
    valuev <- rep("NA", nrow(df.sri)) # default variable values
    cnv[grepl(patti, cnv)]
  })
  names(lmapi) <- names(lgrep)
  lmapi
})
names(lmap) <- names(lsr)
# get mapped values by matched colnames
dfmap <- do.call(rbind, lapply(seq(length(lsr)), function(ii){
  srpi <- names(lsr)[ii]; sri <- lsr[[ii]]
  df.sri <- sri[,grepl("^sra_attribute.*",colnames(sri))]
  cnv <- colnames(df.sri)
  dfmapi <- do.call(cbind, lapply(seq(length(lgrep)), function(ii){
    varname <- names(lgrep)[ii]; patti <- lgrep[[ii]]
    valuev <- rep("NA", nrow(df.sri)) # default variable values
    which.cn <- grepl(patti, cnv)
    # replace valuev if any matching colnames
    if(length(which.cn) > 0){ 
      # uniform formatting of matching values
      valuev <- as.character(apply(df.sri[,which.cn,drop=F], 1, function(ri){
        paste0(gsub("'", "", 
                    gsub(" ", "_", tolower(ri))),
               collapse = ";")
        }))
    }
    dfi <- data.frame(var=valuev); colnames(dfi) <- varname; dfi
  }))
  dfmapi$srp <- srpi; dfmapi$run_id <- sri$external_id; dfmapi
}))

# save mappings
lmd <- list(varname.mappings = lmap, dfmap = dfmap)
lmd.fname <- "list-metadata_dorsolateral-human_recount3.rda"
save(lmd, file = file.path(dir.name, lmd.fname))

#----------------------
# get library prep info
#----------------------
# pattern to search library prep
# get library prep grepl patterns
get_grepl <- function(strv){
  paste0(unlist(lapply(strv, function(ii){
    paste0(".*", ii, ".*")})), collapse = "|")
}
grepl.strv.polya <- c("poly\\(A\\)", "polyA", "poly-A", "poly A")
grepl.strv.ribozero <- c("ribo z", "riboz", "riboZ", "RiboZ")
grepl.strv <- c(grepl.strv.polya, grepl.strv.ribozero)
grepl.pattv <- get_grepl(c(grepl.strv.polya, grepl.strv.ribozero))
grepl.pattv.polya <- get_grepl(grepl.strv.polya)
grepl.pattv.ribozero <- get_grepl(grepl.strv.ribozero)
# table(grepl(grepl.pattv, cd[,libprep.cnv]))
# map library prep info
libprep.cnv <- c("sra.library_construction_protocol", "sra.study_abstract")

df.prep <- do.call(rbind, lapply(seq(length(lcd)), function(ii){
  cdi <- lcd[[ii]]; dfi <- as.data.frame(map_libprep(cdi))
  dfi$project_id <- cdi$study; dfi$run_id <- cdi$external_id
  dfi}))

table(df.prep$libprep.available)
# FALSE  TRUE 
# 250   174

table(df.prep$pa.true)
# FALSE  TRUE 
# 302   122

table(df.prep$rz.true)
# FALSE  TRUE 
# 332    92

table(df.prep$pa.true, df.prep$rz.true)
#       FALSE TRUE
# FALSE   250   52
# TRUE     82   40

#-----------------------
# harmonize library info
#-----------------------
identical(df.prep$run_id, dfmap$run_id) # [1] TRUE
dfall <- cbind(dfmap, df.prep)[,c(1:9,12:16)]

# summarize lib prep variables
table(dfall$library, dfall$pa.true)
#           FALSE TRUE
#             302   82
# polya+       0   20
# ribozero     0   20
table(dfall$library, dfall$rz.true)
#           FALSE TRUE
#             332   52
# polya+       0   20
# ribozero     0   20

# use conditional to get final lib prep info
libvar <- ifelse(dfall$library=="polya"|dfall$pa.true, "polya",
                 ifelse(dfall$library=="ribozero"|dfall$rz.true, 
                        "ribozero", "NA"))
dfall$lib.prep.harmonized <- libvar
table(dfall$lib.prep.harmonized)
# NA    polya ribozero 
# 250      122       52

# save metadata
# dfall
dfall.fname <- "df-md-mapped_lib-prep-harmonized_dorsolateral-human_recount3.rda"
save(dfall, file = file.path(dir.name, dfall.fname))
# lmd, including info
lmd.write <- list(script.fname = "recount3_dataset_properties.R",
                  lmap.colnames = lmap, dfmap = dfmap, dfprep = df.prep, 
                  dfall = dfall)

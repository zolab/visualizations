# ##########################################
# DESCRIPTION: 
# This file imports the targex output file and (i) collapses the results at 
# the note level, (ii) selects a subset of the extracted notes based upon the 
# distribution of notes across time, (iii) visualises and examines the output by
# (a) generating summary statistics, (b) graphing the distribution of notes 
# across time and patients, (b) graphing the distribution of notes across time 
# for a subset of individual patients, (c) graphing the distbution of concept 
# mentions across time dor a susbet of indiviual patients and paring that 
# information with encounter, ed encounter, diagnoses and procedure information 
# (to validate the results), (iv) formats the targex output file according to 
# the specifications of the colouring macro
# ##########################################
# Creator: Clara Marquardt
# Date: 1st December
# ##########################################
# Language: R
# ##########################################
# TO-DO-LIST 
# FORMAT
## [] comment code
## [] graphs - improve formatting & ensure that formatting of e.g. geom_tiles 
## work/look acceptable for any number of concepts 
## multiple geom_tiles on one page - manually control scale gradient cut-offs  
## -- COLOURS ARE TOO BRIGHT (I.E LIGHT)
### so that standarized across graphs
### [] correlation matrix - add lines to separate concepts
## [] save all graphs as separate png files alongside the pdf output
# GENERALISE
## [] working directory - currently nothing is saved as an intermediary ouput
## [] ensure that raw files also work with non-pants (i.e. non bwh_ed) versions 
## of these files (problem only - if at all - for ed encoutner file - all other 
## files = raw files)
### [] output - currently saved as csv file vs. also explored option of using 
### tableGrob(table name) combined wiht gridarrange(plot1, table, asTable=T)
## [] implement call to VB macro - if not from within R then from within terminal - 
## i.e. save call to terminal at the end of the code file
## [] loading the full input file - assuming that called notes.chunks.sent in 
## workspace & that 157 chunks
## [] ensure that load formatted diagnoses and procedure files - rather than 
## formatted but already subset procedure files & ensure that current forattign 
## is correct
# QUESTIONS
## [] time frame across which to map/examine & how to specify this dynamically & 
### how to deal with leap years
### [] types of data sources that are relevant enough to be mapped across projects 
### (rather than being specific to the pants project) incl. help files such as 
### diagnoses/procedure codes
### [] combine the subsetting and visualisation file in one or split into two 
### chunks? Subsetting - does the current approach make sense - i.e. basing it 
### on the length of the period and the length between notes? 
### [] sensibility of validation by other data files - NLP as subset vs. superset
### [] currently sum concepts across sentences and sum concept member dummies 
### (converting them into count variables rather than dummies) - then sum both 
### across notes to  generate weekly counts - problem: concept member dummies 
### generated as dummies not count variables
### [] output - currently saved as csv file vs. also explored option of using 
### tableGrob(table name) combined wiht gridarrange(plot1, table, asTable=T)
### [] subset by number of concept member drivers - how to best implement this & 
### alter implementation once all concept members are laoded
## [] hierarchy of restrictins - (1) stability (2) activity (which in turn (a) 
## gap betwen notes and (b) duration of note periods) -- if this implementation is 
## to be used fix code so htat implement int his order while retaining them 
## seperateld y (e.f. active but nonstable) - note current code uses the 
## combination of stability and activity when defininf active notes
################################################################################


################################################################################
######################  (A) SET-UP - 'Control File'  ###########################
################################################################################
### operational specifications
# (a) working directory
setwd("/data/zolab/pants") 
modified_data_folder <- "modified_data/"
output_data_folder <- "output/"
filename <- "_cancerchemo"
currentdate <- as.character(format(Sys.time(), "%d-%m-%Y")) 
# (b) libraries and functions
source("/data/zolab/pants/helpers/30_day_cancer_chemo/pants_cancer_chemo_functions.R")
source("/data/zolab/pants/helpers/30_day_cancer_chemo/pants_cancer_chemo_libraries.R")

# (c) targex output file
output_filename <- 
  "/data/zolab/pants/modified_data/cancerchemo_chemodrug_targex_output_mod_cancerchemo.csv"

# (d) other raw data sources - load diagnosis/procedure files that are formatted
# TO_DO: Figure out how to load raw data then format in code - for now load formatted 
# and subsetted file and comemnt out subsetting 
encounter_filename <- "/data/zolab/pants/raw_data/rpdr_2010_2012_Enc.csv"
ed_encounter_filename <- "/data/zolab/pants/modified_data/ed_encounters_mod.csv"
diagnoses_filename <- "/data/zolab/pants/raw_data/rpdr_2010_2012_Dia.csv"
procedures_filename <- "/data/zolab/pants/raw_data/rpdr_2010_2012_Prc.csv"

# (e) full input (e.g. full notes) file (chunkified - notes only)
full_input_file_name <- 
  "/data/zolab/bwh_ed/targex/lno/store/bwh100k_lno_chunk_all_sent.RData"
# NOTE: Assumigng that when loaded the file will be called "notes.chunks.sent" & 
# containing 157 chunks

################################################################################
### inputs

# (0) number of concept members that need to be included for a note to be even 
# considered (pre all other analysis) (if no subsetting based on this - set to NA)
min_concept_members <- 2

# (a) variables to determine which (if not all) notes are selected from among the 
# entire set of extracted notes 
time_frame_frequency_weeks <- 12  # max gap (weeks) between notes & extension of 
							      # 'period' after last note in a sequence
time_frame_period_weeks <- 12  # min length (weeks) over which a sequence of 
							   # contigious notes has to extend

# XXX Note: 3 options - (a) all notes - set both variables = NA, (b) notes 
# separated by max xxx weeks regardless of the overall duration -  set 
# time_frame_frequency_weeks= number & time_frame_period_weeks = 0, (c) notes 
# separated by xxx weeks and extending over at least xxx weeks - set both 
# variables = number
# XXX Note: If subsetting is implementing - encounters, ed encounters, diagnoses 
# and procedures are deemed 'relevant' if they fall anywehre within a sequence 
# of notes (date of first note in a sequence to date of last note in a sequence + 
# 'time_frame_frequency_weeks' weeks)
# XXX Note: All time periods between notes (and encounters, etc.) are determined 
# based on the distance between the weeks of eachn otes, i.e. distance between 
# mondays of hte weeks in which notes fall

# (b) variables to determine the type of encounter (clinic cat) that will be 
# examined
enc_clinic_cat <- "onc"
# XXX Note: See clinic x_walk for the list of possible clinic types - 
# to examine all encounters - input NA

# (c) variables to determine the types of diagnoses and procedures that will be 
# examined
procedure_codes_filename <- 
  "/data/zolab/pants/helpers/30_day_cancer_chemo/chemo_codes_mod.csv" 
diagnosis_codes_filename <-
  "/data/zolab/pants/helpers/30_day_cancer_chemo/cancer_codes_mod.csv" 

# (d) variables to determine whether concepts or concepts and concept members 
# are mapped 
variables_to_graph<- "concept_list_expanded"
# XXX Note: Input 'concept_list' or 'concept_list_expanded'


################################################################################
######################  (A) SET-UP - Load Files ################################
################################################################################
################################################################################
# (a) load targex output file & rename (for convenience)
output<- setkey(fread(output_filename), empi)
setnames(output, c("noteid", "sentid", "senttext"), c("note_id", "sent_id", 
  "sent_text"))


# (b) load raw data & specify column classes & format dates/diagnoses & 
# procedure codes

# encounters
encounters <- setkey(fread(encounter_filename, colClasses=c(diagnosis_7="character", 
  diagnosis_8="character",diagnosis_9="character",diagnosis_10="character", 
  admit_source="character", empi="character")), empi)

# ed encounters
ed_encounters <- setkey(fread(ed_encounter_filename,colClasses=c(empi="character")), 
  empi)

# diagnoses  (pull in file which has standarised codes)
#diagnoses_raw<- setkey(fread(diagnoses_filename, colClasses=c(code="character",
#  empi="character")), empi)
diagnoses_raw<- setkey(fread("/data/zolab/pants/modified_data/diagnoses_cancer.csv", 
  colClasses=c(code="character", empi="character")), empi)

# procedures (pull in file which has standarised codes)
#procedures_raw <- setkey(fread(procedures_filename,colClasses=c(code="character",
#  procedure_flag="character", empi="character")),empi)
procedures_raw<- setkey(fread("/data/zolab/pants/modified_data/procedures_chemo.csv", 
  colClasses=c(code="character", empi="character")), empi)

# (c) load help files
clinic_xwalk <- fread("/data/zolab/pants/helpers/rpdr_clinic_xwalk_15_07_14.csv", 
  colClasses=c(clinic_no="character"))

# pull in code files which are standarised
procedure_codes <- fread(procedure_codes_filename)

diagnoses_codes <- fread(diagnosis_codes_filename)

# (d) load the full notes (input) - this is only necessary for the colouring macro 
# output
load(file = full_input_file_name)

# (d) load the full notes (strucutred) & assign note id 
notes_data_structured <- fread(full_input_structured_name)[, note_id:=1:.N]

################################################################################
############  (B) FORMAT THE TARGEX OUTPUT FILE (COLLAPSE & SUBSET)  ###########
################################################################################
### collapse the targex output file at the note level (for each concept - sentence 
### and dummy variable for each note - to be used as input into colouring macro)

# (a) extract the concept_list/concept_member list
concept_list <- copy(names(output)[!(names(output) %like% 
  "\\.| ")])
concept_list_expanded <- copy(names(output)[names(output) 
  %like% "[A-Z]"])
concept_member_list<- setdiff(concept_list_expanded,concept_list)


# (b) gather all concept variables into two columns and drop the 
# original concept variables
output <- as.data.table(gather(output, concept_name, concept_count, -(1:5), 
	-matches("\\.")))

# set the sentence variable = "" for rows in which the concept count variable=0
output[concept_count==0, sent_text:=""]

# generate the concept_dummy and concept_sentence variables
output <- output%>%
  unite(temp, concept_count, sent_text,sep="------") %>%
  spread(concept_name, temp)
output <- as.data.table(output)

for (var in concept_list) {
  var_name <- parse(text=var)
  output[,c(paste(rep(var,2),c("dummy","sentence"),sep="_")):=
  tstrsplit(eval(var_name), "------")]
}

output[, c(concept_list):=NULL]

# prepare to collapse - sum all concept dummies, collapse all sentences and sum 
# concept member dummies (retain dummy format) (note: format of sentences is 
# determined by requirements of colouring macro)
output[, c(paste0(concept_list, "_dummy_sum")) := 
  lapply(.SD, function(x) sum(as.numeric(x))), by=c("note_id"), 
  .SDcols=paste0(concept_list,"_dummy")][, c(paste0(concept_list, "_dummy")):=
  NULL]

output[, c(paste0(concept_member_list,"_dummy_sum")) := 
  lapply(.SD, function(x) sum(as.numeric(x))), by=c("note_id"), 
  .SDcols=concept_member_list][, c(concept_member_list):=
  NULL]

# alternative - retain dummy format 
'output[, c(paste0(concept_member_list,"_dummy_sum")) := 
  lapply(.SD, function(x) ifelse(sum(as.numeric(x))>0,1,0)), by=c("note_id"), 
  .SDcols=concept_member_list][, c(concept_member_list):=
  NULL]'
# 

for (var in c(paste(concept_list, "sentence",sep="_"))) {
  var_name <- parse(text=var)
  output[is.na(eval(var_name)),  c(paste(var)) :=.("")]
  output[,  c(paste(var)) :=
  paste0("...", paste(eval(var_name),collapse="..."), "..."), by=c("note_id")]
}

# collapse by note_id
output <- unique(output, by=c("note_id"))

################################################################################
### generate variables that will be employed in subsequent analysis steps & subset 
### set of relevant notes based on specificaion of across time requirements

# list of 'core' variables
variables <- copy(names(output))

# time windows, year, month, week  -- week number (%W) & day number (%u) -  
# first monday in a month = W1 & monday = day 1
output[, c("lno_year_week") := .(paste0(format(strptime(lno_date, "%Y-%m-%d"), 
  "%Y-%W"), "-1"))]
output[, c("lno_year", "lno_month", "lno_week") := 
  .(format(strptime(lno_date, "%Y-%m-%d"), "%Y"),format(strptime(lno_date, 
  "%Y-%m-%d"), "%m"), format(strptime(lno_date, "%Y-%m-%d"), "%W"))]

# XXX TO-DO: Resolve this - how to deal with leap years & week 0 observations
# move all week 0 and week 53 observations into week 1 and week 52 
output[lno_week=="00", c("lno_week","lno_year_week"):=.("01", gsub("-00-", "-01-", 
  lno_year_week))]
output[lno_week=="53", c("lno_week","lno_year_week"):=.("52", gsub("-53-", "-52-", 
  lno_year_week))]

# identify if concept is uniquely driving the identification of a note 
# unique driver  if concept count = 1
output[, c(paste0("temp_", 1:length(concept_list))) :=lapply(.SD, 
  function(x) x!=0), .SDcols=c(paste0(concept_list, "_dummy_sum"))]
output[, concept_count:=.(rowSums(.SD)), .SDcols=c(paste0("temp_", 
  1:length(concept_list)))][, c(paste0("temp_", 1:length(concept_list))):= NULL]

output[, c(paste0("temp_", 1:length(concept_list_expanded))) :=lapply(.SD, 
  function(x) x!=0), .SDcols=c(paste0(concept_list_expanded, "_dummy_sum"))]
output[, concept_count_mod:=.(rowSums(.SD)), .SDcols=c(paste0("temp_", 
  1:length(concept_list_expanded)))][, c(paste0("temp_", 1:length(concept_list_expanded))):= NULL]


# identify if concept member is uniquely driving the identification of a note 
# susbet based on this if this is to be implemented
# XXX TO_DO: Decide how to best implement this depend on file structre - how granular a 
# division into concept members & change above to concept member list 
# rather than concept list
#output[, c(paste0("temp_", 1:length(concept_member_list))) :=lapply(.SD, 
#  function(x) x!=0), .SDcols=c(paste0(concept_member_list, "_dummy_sum"))]
#output[, concept_member_count:=.(rowSums(.SD)), .SDcols=c(paste0("temp_", 
#  1:length(concept_member_list)))][, c(paste0("temp_", 
#  1:length(concept_member_list))):= NULL]
output[, concept_member_count:=.(rowSums(.SD)), .SDcols=c(
  paste0(concept_member_list, "_dummy_sum"))]

if (!is.na(min_concept_members)){
  output[concept_member_count<min_concept_members, stable_identification:=0][
    is.na(stable_identification), stable_identification:=1]
} else {
  output[, stable_identification:=1]
}

# determine the number of notes per patient
output[, notes_per_patient := .N, by=c("empi")]

# check if subsetting is to be implemented & if yes - implement
if (!is.na(time_frame_frequency_weeks)) {

  output <- output[order(empi, lno_date)]
  # 1 note/patient - nonactive
  output <- output[notes_per_patient==1, active_note_broad:=0] 
  # non-stable - non-active 
  output <- output[stable_identification==0, active_note_broad:=0] 
  output[notes_per_patient!=1 & stable_identification!=0, 
    days_from_last_note_lag:=round(as.numeric((
    strptime(lno_year_week, "%Y-%W-%u")-strptime(shift(lno_year_week, 1L, 
    type="lag"), "%Y-%W-%u")), units = "weeks"),digits=0),by=c("empi")]
  output[notes_per_patient!=1 & stable_identification!=0, 
    days_from_last_note_lead:=round(as.numeric((
    strptime(lno_year_week, "%Y-%W-%u")-strptime(shift(lno_year_week, 1L, 
    type="lead"), "%Y-%W-%u")), units = "weeks"),digits=0),by=c("empi")]

  # merge time to previous and time to next nore and identify active notes
  output[, days_from_last_note:=days_from_last_note_lag][
    is.na(days_from_last_note_lag) & is.na(active_note_broad), days_from_last_note:=
    days_from_last_note_lead][days_from_last_note_lag > 
    abs(days_from_last_note_lead), days_from_last_note:=days_from_last_note_lead]

  output[abs(days_from_last_note) > time_frame_frequency_weeks, active_note_broad:=0][
    is.na(active_note_broad), active_note_broad:=1]

  # generate identifiers to identify periods of active onc status
  # last observation
  output[is.na(days_from_last_note_lead) & active_note_broad==1, note_period:=1] 
  # 1st observation
  output[is.na(days_from_last_note_lag) & active_note_broad==1 , note_period:=0] 
  # end of period
  output[abs(days_from_last_note_lead) > time_frame_frequency_weeks & 
    active_note_broad==1 & is.na(note_period), note_period:=2] 
  # beginning of period
  output[days_from_last_note_lag > time_frame_frequency_weeks & 
    active_note_broad==1 & is.na(note_period), note_period:=3] 

  # gourp into periods and assign note_period index 
  output[note_period %in% c(0,3), note_period_index:=1:.N]
  output[active_note_broad==1,note_period_index:=na.locf(note_period_index)]

  # determine number of unique dates wihtn each note period 
  output[active_note_broad==1, note_period_week_count:=length(unique(lno_year_week)), 
    by=c("note_period_index")]

  # determine length of each note period & calculate the extended note period end
  output[active_note_broad==1 , note_period_days_end:=.SD[c(.N)], 
    by=c("note_period_index"), .SDcols=c("lno_year_week")]
  output[active_note_broad==1 , note_period_days_beg:=.SD[c(1)], 
    by=c("note_period_index"), .SDcols=c("lno_year_week")]
  output[active_note_broad==1 , note_period_days:=
    round(as.numeric((strptime(note_period_days_end, "%Y-%W-%u")-
    strptime(note_period_days_beg, "%Y-%W-%u")), units = "weeks"),digits=0)]
  output[, note_period_days_end_extended:=format(as.IDate(
    strptime(note_period_days_end, "%Y-%W-%u")+weeks(1)), "%Y-%W-%u")]

  # drop variables
  output[, c("days_from_last_note_lag", "days_from_last_note_lead"):=NULL]

  # refine active notes 
  output[, active_note:=active_note_broad][note_period_days < 
    time_frame_period_weeks, active_note:=0]

  # determine hte active notes per patient 
  output[active_note!=0, active_notes_per_patient := .N, by=c("empi")][is.na(
  active_notes_per_patient), active_notes_per_patient:=0]

} else {
  # if subsetting is not to be implemented - manually deem all notes active
  output[, active_note:=1]
  output[, active_note_broad:=1]
  output[, active_notes_per_patient:=notes_per_patient]
}


# define the lsit of additional variables
additional_variables <- setdiff(names(output), variables)

# SAVE - DATA OUTPUT 1 - Targex output file collapsed at the note level with 
# added variables & subsettign
write.csv(output, file=paste0(modified_data_folder, "formatted_targex_output", 
  filename,".csv"), row.names=F)


################################################################################
############  (C) GENERATE PATIENT SUBSETS (FOR THE VISUAL ANALYSIS)  ##########
################################################################################
set.seed(0)

# (a) 50 patient subset - balanced in terms of the active number of notes per 
# patient 
empi_index_expanded <- unique(output[,.(empi, active_notes_per_patient)], 
  by=c("empi"))
subset <- with(empi_index_expanded, createDataPartition(active_notes_per_patient, 
  p=50/nrow(empi_index_expanded)))
empi_index_expanded <- empi_index_expanded[subset$Resample1]
empi_index_expanded <- empi_index_expanded[1:50] 

# 300 patient subset for colour coding macro - balanced across number of 
# active notes but only people with a least one active note
empi_index_expanded_large <- unique(output[!active_notes_per_patient==0,
  .(empi, active_notes_per_patient)], by=c("empi"))
subset <- with(empi_index_expanded_large, createDataPartition(
  active_notes_per_patient, p=300/nrow(empi_index_expanded_large)))
empi_index_expanded_large <- empi_index_expanded_large[subset$Resample1]
empi_index_expanded_large <- empi_index_expanded_large[1:300] 

# 20 patient subset -  balanced in terms of the number of active notes per 
# patient 
empi_index_small <- unique(output[active_notes_per_patient!=0, 
  .(empi, active_notes_per_patient)])
subset <- with(empi_index_small, createDataPartition(active_notes_per_patient, 
  p=20/nrow(empi_index_small)))
empi_index_small <- empi_index_small[subset$Resample1]
empi_index_small <- empi_index_small[1:20] 

# 20 patient subset - 20 patients with the highest number of active notes per 
# patient
empi_index_high_frequency <- unique(output[active_notes_per_patient!=0, 
  .(empi, active_notes_per_patient)][order(-active_notes_per_patient)])
empi_index_high_frequency <- empi_index_high_frequency[1:20] 

################################################################################
################  (D) FORMAT OTHER FILES (ENC/ED ENC/DIAG/PROC)  ###############
################################################################################
### encounters

# merge with rpdr clinic crosswalk 
encounters <- encounters[clinic_xwalk, .(empi, clinic_cat, hospital, admit_date),
   on=c("clinic_name"), nomatch=NA][, encounter_number:=1:.N]

# subset to specified clinic type & patients with identified note
if (!is.na(enc_clinic_cat)) {
  output_encounters <- encounters[clinic_cat==enc_clinic_cat & 
    empi %in% output$empi]
} else {
  output_encounters <- encounters[empi %in% output$empi]
}

# generate year - week groups 
output_encounters[, c("admit_year_week") := .(paste0(format(strptime(admit_date, 
  "%m/%d/%Y"), "%Y-%W"), "-1"))]
output_encounters[, c("admit_year", "admit_month", "admit_week") := 
  .(format(strptime(admit_date, "%m/%d/%Y"), "%Y"),
  format(strptime(admit_date, "%m/%d/%Y"), "%m"), format(strptime(admit_date, 
  "%m/%d/%Y"), "%W"))]

# move all week 0 and week 53 observations into week 1 and week 52 
output_encounters[admit_week=="00", c("admit_week","admit_year_week"):=
  .("01", gsub("-00-", "-01-", admit_year_week))]
output_encounters[admit_week=="53", c("admit_week","admit_year_week"):=
  .("52", gsub("-53-", "-52-", admit_year_week))]


# if subsetting is to be implemented mark all encounters asa relevant which fall 
# within note periods - otherwise mark all encoutners as relevant
if (!is.na(time_frame_frequency_weeks)) {
  output_encounters_exp <- output_encounters[output, .(empi, admit_date, 
    admit_year_week, encounter_number, note_id, lno_date, lno_year_week, 
    note_period_index, note_period_days_beg, note_period_days_end, 
    note_period_days_end_extended, active_note), on=c("empi"), nomatch=NA, 
    allow.cartesian=T]
  output_encounters_exp[active_note==1 & strptime(admit_year_week, "%Y-%W-%u")>=
    strptime(note_period_days_beg, "%Y-%W-%u") & strptime(admit_year_week, 
    "%Y-%W-%u") <= strptime(note_period_days_end_extended, "%Y-%W-%u"), 
    relevant_enc:=1][is.na(relevant_enc),relevant_enc:=0]
  output_encounters <- unique(output_encounters_exp[, .(empi, admit_date, 
    admit_year_week, encounter_number, relevant_enc)][order(encounter_number, 
    -relevant_enc)], by=c("encounter_number"))
} else {
  # mark all encounters as relevant
  output_encounters[, relevant_enc:=1]
}

output_encounters_coll <- output_encounters[,':='(enc_count=.N, 
  enc_count_active=sum(relevant_enc==1)), by=c("empi", "admit_year_week")]
output_encounters_coll <- as.data.table(rbind(output_encounters_coll[, 
  .(empi, admit_year_week, enc_count, cat="all")], setnames(
  output_encounters_coll[, .(empi, admit_year_week, enc_count_active, 
  cat="active")], names(output_encounters_coll[, .(empi, admit_year_week, enc_count, 
  cat="all")]))))

# SAVE - DATA OUTPUT 2 - Encounters for patients with identified notes linked to 
# these ntoes by date
write.csv(output_encounters_coll, file=paste0(modified_data_folder,
  "targex_output_encounters", filename,".csv"), row.names=F)

################################################################################
### ED encounters

# susbet to encounters with identified notes
output_ed_encounters <- ed_encounters[empi %in% output$empi][, .(empi, ed_date, 
  encounter_number)]

# generate year - week groups 
output_ed_encounters[, c("ed_year_week") := paste0(format(strptime
  (ed_date, "%d%b%Y"), "%Y-%W"), "-1")]
output_ed_encounters[, c("ed_year", "ed_month", "ed_week") := 
  .(format(strptime(ed_date, "%d%b%Y"), "%Y"),
  format(strptime(ed_date, "%d%b%Y"), "%m"), format(strptime(ed_date, 
  "%d%b%Y"), "%W"))]

# move all week 0 and week 53 observations into week 1 and week 52 
output_ed_encounters[ed_week=="00", c("ed_week","ed_year_week"):=
  .("01", gsub("-00-", "-01-", ed_year_week))]
output_ed_encounters[ed_week=="53", c("ed_week","ed_year_week"):=
  .("52", gsub("-53-", "-52-", ed_year_week))]

# if subsetting is to be implemented mark all encounters asa relevant which fall 
# within note periods - otherwise mark all encoutners as relevant
if (!is.na(time_frame_frequency_weeks)) {
  output_ed_encounters_exp <- output_ed_encounters[output, .(empi, ed_date, 
    ed_year_week, encounter_number, note_id, lno_date, lno_year_week, note_period_index, 
    note_period_days_beg, note_period_days_end, note_period_days_end_extended, 
    active_note), on=c("empi"), nomatch=NA, allow.cartesian=T]
  output_ed_encounters_exp[active_note==1 & strptime(ed_year_week, "%Y-%W-%u")>=
    strptime(note_period_days_beg, "%Y-%W-%u") & strptime(ed_year_week, 
    "%Y-%W-%u") <= strptime(note_period_days_end_extended, "%Y-%W-%u"),
    relevant_enc:=1][is.na(relevant_enc),relevant_enc:=0]
  output_ed_encounters <- unique(output_ed_encounters_exp[, .(empi, ed_date, 
    ed_year_week,encounter_number, relevant_enc)][order(encounter_number, 
    -relevant_enc)], by=c("encounter_number"))
} else {
  # mark all encounters as relevant
  output_ed_encounters[, relevant_enc:=1]
}

output_ed_encounters_coll<- output_ed_encounters[,':='(ed_enc_count=.N, 
  ed_enc_count_active=sum(relevant_enc==1)), by=c("empi", "ed_year_week")]
output_ed_encounters_coll <- as.data.table(rbind(output_ed_encounters[, 
  .(empi, ed_year_week, ed_enc_count, cat="all")], setnames(
  output_ed_encounters[, .(empi, ed_year_week, ed_enc_count_active, 
  cat="active")], names(output_ed_encounters[, .(empi, ed_year_week, ed_enc_count, 
  cat="all")]))))

# SAVE - DATA OUTPUT 3 - ED encounters for patients with identified notes linked to 
# these ntoes by date
write.csv(output_ed_encounters_coll, file=paste0(modified_data_folder,
  "targex_output_ed_encounters", filename,".csv"), row.names=F)


################################################################################
### Diagnoses

# extract the relevant diagnoses & for the relevant patients
# output_diagnoses <- diagnoses_raw[diagnoses_codes, on=c(code_type="code_type",
# code="cancer_code"), .(empi, mrn, mrn_type, diagnosis_name, code_type, code,
# icd_date,clinic,hospital), nomatch=0][empi %in% output$empi][, 
#  encounter_number:=1:.N]
output_diagnoses <- diagnoses_raw[, .(empi, mrn, mrn_type, diagnosis_name, 
  code_type, code, icd_date,clinic,hospital)][empi %in% output$empi][, 
encounter_number:=1:.N]

# generate year - week groups 
output_diagnoses[, c("diagnosis_year_week") := paste0(format(strptime
  (icd_date, "%d%b%Y"), "%Y-%W"), "-1")]
output_diagnoses[, c("diagnosis_year", "diagnosis_month", "diagnosis_week") := 
  .(format(strptime(icd_date, "%d%b%Y"), "%Y"),
  format(strptime(icd_date, "%d%b%Y"), "%m"), format(strptime(icd_date, 
  "%d%b%Y"), "%W"))]

# move all week 0 and week 53 observations into week 1 and week 52 
output_diagnoses[diagnosis_week=="00", c("diagnosis_week","diagnosis_year_week"):=
  .("01", gsub("-00-", "-01-", diagnosis_year_week))]
output_diagnoses[diagnosis_week=="53", c("diagnosis_week","diagnosis_year_week"):=
  .("52", gsub("-53-", "-52-", diagnosis_year_week))]

# if subsetting is to be implemented mark all diagnoses asa relevant which fall 
# within note periods - otherwise mark all diagnsoes as relevant
if (!is.na(time_frame_frequency_weeks)) {
  output_diagnoses_exp <- output_diagnoses[output, .(empi, icd_date, 
    diagnosis_year_week, encounter_number, note_id, lno_date, lno_year_week, 
    note_period_index, note_period_days_beg, note_period_days_end, 
    note_period_days_end_extended, 
    active_note), on=c("empi"), nomatch=NA, allow.cartesian=T]
  output_diagnoses_exp[active_note==1 & strptime(diagnosis_year_week, "%Y-%W-%u")>=
    strptime(note_period_days_beg, "%Y-%W-%u") & strptime(diagnosis_year_week, 
    "%Y-%W-%u") <= strptime(note_period_days_end_extended, "%Y-%W-%u"),
    relevant_enc:=1][is.na(relevant_enc),relevant_enc:=0]
  output_diagnoses <- unique(output_diagnoses_exp[, .(empi, icd_date, 
    diagnosis_year_week,encounter_number, relevant_enc)][order(encounter_number, 
    -relevant_enc)], by=c("encounter_number"))
} else {
  # mark all diagnoses as relevant
  output_diagnoses[, relevant_enc:=1]
}

output_diagnoses_coll<- output_diagnoses[,':='(diagnosis_count=.N, 
  diagnosis_count_active=sum(relevant_enc==1)), by=c("empi", "diagnosis_year_week")]
output_diagnoses_coll <- as.data.table(rbind(output_diagnoses[, 
  .(empi, diagnosis_year_week, diagnosis_count, cat="all")], setnames(
  output_diagnoses[, .(empi, diagnosis_year_week, diagnosis_count_active, 
  cat="active")], names(output_diagnoses[, .(empi, diagnosis_year_week, 
  diagnosis_count, cat="all")]))))

# SAVE - DATA OUTPUT 4 - Diagnsoes for patients with identified notes linked to 
# these ntoes by date
write.csv(output_diagnoses_coll, file=paste0(modified_data_folder,
  "targex_output_diagnoses", filename,".csv"), row.names=F)


################################################################################
### Procedures

# extract the relevant procedures & for the relevant patients
#output_procedures<- procedures_raw[procedure_codes, on=c(code_type="code_type",
#  code="chemo_code"), .(empi, mrn, mrn_type, procedure_name, code_type, code, date,
#  clinic,hospital),nomatch=0][empi %in% output$empi][, 
#  encounter_number:=1:.N]
output_procedures<- procedures_raw[, .(empi, mrn, mrn_type, procedure_name, 
  code_type, code, date,clinic,hospital)][empi %in% output$empi][, 
  encounter_number:=1:.N]

# generate year - week groups 
output_procedures[, c("procedure_year_week") := paste0(format(strptime
  (date, "%d%b%Y"), "%Y-%W"), "-1")]
output_procedures[, c("procedure_year", "procedure_month", "procedure_week") := 
  .(format(strptime(date, "%d%b%Y"), "%Y"),
  format(strptime(date, "%d%b%Y"), "%m"), format(strptime(date, 
  "%d%b%Y"), "%W"))]

# move all week 0 and week 53 observations into week 1 and week 52 
output_procedures[procedure_week=="00", c("procedure_week","procedure_year_week"):=
  .("01", gsub("-00-", "-01-", procedure_year_week))]
output_procedures[procedure_week=="53", c("procedure_week","procedure_year_week"):=
  .("52", gsub("-53-", "-52-", procedure_year_week))]

# if subsetting is to be implemented mark all procedures asa relevant which fall 
# within note periods - otherwise mark all procedures as relevant
if (!is.na(time_frame_frequency_weeks)) {
  output_procedures_exp <- output_procedures[output, .(empi, date, 
    procedure_year_week, encounter_number, note_id, lno_date, lno_year_week, 
    note_period_index, note_period_days_beg, note_period_days_end, 
    note_period_days_end_extended, active_note), on=c("empi"), nomatch=NA, 
    allow.cartesian=T]
  output_procedures_exp[active_note==1 & strptime(procedure_year_week, "%Y-%W-%u")>=
    strptime(note_period_days_beg, "%Y-%W-%u") & strptime(procedure_year_week, 
    "%Y-%W-%u") <= strptime(note_period_days_end_extended, "%Y-%W-%u"),
    relevant_enc:=1][is.na(relevant_enc),relevant_enc:=0]
  output_procedures <- unique(output_procedures_exp[, .(empi, date, 
    procedure_year_week,encounter_number, relevant_enc)][order(encounter_number, 
    -relevant_enc)], by=c("encounter_number"))
} else {
  # mark all procedures as relevant
  output_procedures[, relevant_enc:=1]
}

output_procedures_coll<- output_procedures[,':='(procedure_count=.N, 
  procedure_count_active=sum(relevant_enc==1)), by=c("empi", "procedure_year_week")]
output_procedures_coll <- as.data.table(rbind(output_procedures[, 
  .(empi, procedure_year_week, procedure_count, cat="all")], setnames(
  output_procedures[, .(empi, procedure_year_week, procedure_count_active, 
  cat="active")], names(output_procedures[, .(empi, procedure_year_week, 
  procedure_count, cat="all")]))))

# SAVE - DATA OUTPUT 5 - Procedures for patients with identified notes linked to 
# these ntoes by date
write.csv(output_procedures_coll, file=paste0(modified_data_folder,
  "targex_output_procedures", filename,".csv"), row.names=F)


################################################################################
####################### (E) OUTPUTS  I  - SUMMARY STATS ########################
################################################################################
### Summary Stats

# main summary stars
targex_output_summary_stats <- list(
  "unique_patients_count" = length(unique(output$empi)),
  "unique_notes_count" = length(unique(output$note_id)), 
  "unique_note_dates" = length(unique(output$lno_date)), 
  "unique_note_patient_dates" = nrow(unique(output[,.(lno_date, empi)],by=NULL)), 
  "unique_onc_periods" = NA
) 

# active ntote specific summary stats if subsetting is implemented
if (!is.na(time_frame_frequency_weeks)) {

targex_output_summary_stats_stable <- list(
  "unique_patients_count" = length(unique(output[stable_identification==1]$empi)),
  "unique_notes_count" = length(unique(output[stable_identification==1]$note_id)), 
  "unique_note_dates" = length(unique(output[stable_identification==1]$lno_date)),
  "unique_note_patient_dates" = nrow(unique(output[stable_identification==1,.(lno_date, empi)],by=NULL)), 
  "unique_onc_periods" = NA
) 

#targex_output_summary_stats_active_nonstable <- list(
 # "unique_patients_count" = length(unique(output[active_note_nonstable!=0]$empi)),
  #"unique_notes_count" = length(unique(output[active_note_nonstable!=0]$note_id)), 
  #"unique_note_dates" = length(unique(output[active_note_nonstable!=0]$lno_date)),
  #"unique_onc_periods" = NA
#) 

targex_output_summary_stats_active <- list(
  "unique_patients_count" = length(unique(output[active_note!=0]$empi)),
  "unique_notes_count" = length(unique(output[active_note!=0]$note_id)), 
  "unique_note_dates" = length(unique(output[active_note!=0]$lno_date)),
  "unique_note_patient_dates" = nrow(unique(output[active_note!=0,.(lno_date, empi)],by=NULL)), 
  "unique_onc_periods" = length(unique(output[active_note==1]$note_period_index))
) 


targex_output_summary_stats<-as.data.frame(t(rbind(as.data.frame(
  targex_output_summary_stats), as.data.frame(
  targex_output_summary_stats_stable),
  as.data.frame(
  targex_output_summary_stats_active))))
targex_output_summary_stats$V0 <- rownames(targex_output_summary_stats)
colnames(targex_output_summary_stats) <- c("all extracted notes", 
  "stable extracted notes", 
  "stable and active extracted notes", "var_names")
targex_output_summary_stats <- as.data.table(targex_output_summary_stats)

} else {
  targex_output_summary_stats<-as.data.frame(t(
    as.data.frame(targex_output_summary_stats)))
  targex_output_summary_stats$V0 <- rownames(targex_output_summary_stats)
  colnames(targex_output_summary_stats) <- c("var", "var_name")
  targex_output_summary_stats <- as.data.table(targex_output_summary_stats)
}

# SAVE - FINAL OUTPUT 1a - Summary stats at the note level
write.csv(targex_output_summary_stats, file=paste0(output_data_folder,
  "targex_output_summary_stats", filename,".csv"), row.names=F)


# Summary stats based on linking notes with encoutners, etc. - only if subsetting 
# is implemented

if (!is.na(time_frame_frequency_weeks)) {
  targex_output_summary_stats_validation <- list(
    "number of patients with a diagnosis/procedure/encounter at any time" = 
     length(unique(c(output_procedures_exp[!is.na(encounter_number) & 
      active_note==1]$empi,output_diagnoses_exp[!is.na(encounter_number)
      & active_note==1]$empi, output_encounters_exp[!is.na(encounter_number) & 
      active_note==1]$empi))), 
    "percentage of patients with a diagnosis/procedure/encounter at any time" = round((
     length(unique(c(output_procedures_exp[!is.na(encounter_number) & 
      active_note==1]$empi,output_diagnoses_exp[!is.na(encounter_number)
      & active_note==1]$empi, output_encounters_exp[!is.na(encounter_number) & 
      active_note==1]$empi)))/length(unique(output[active_note!=0]$empi)))*100, 
      digits=0),

   "percentage of patients with a diagnosis at any time" = round((
     length(unique(output_diagnoses_exp[!is.na(encounter_number)
      & active_note==1]$empi))/length(unique(output[active_note!=0]$empi)))*100, 
      digits=0),
    "percentage of patients with a procedure at any time" = round((
     length(unique(output_procedures_exp[!is.na(encounter_number) & 
      active_note==1]$empi))/length(unique(output[active_note!=0]$empi)))*100, 
      digits=0),
    "percentage of patients with an encounter at any time" = round((
     length(unique(output_encounters_exp[!is.na(encounter_number) & 
      active_note==1]$empi))/length(unique(output[active_note!=0]$empi)))*100, 
      digits=0),
  
    "number of note periods validated by an encounter" = 
      length(unique(output_encounters_exp[relevant_enc==1]$note_period_index)),
    "percentage of note periods validated by an encounter" = round((
       length(unique(output_encounters_exp[relevant_enc==1]$note_period_index))/
       length(unique(output[active_note==1]$note_period_index)))*100, digits=0),
    "number of notes validated by an encounter" = 
       length(unique(output_encounters_exp[relevant_enc==1]$note_id)),     
    "percentage of notes validated by an encounter" = round((
       length(unique(output_encounters_exp[relevant_enc==1]$note_id))/
       length(unique(output[active_note==1]$note_id)))*100, digits=0),

    "number of note periods validated by a procedure" = 
      length(unique(output_procedures_exp[relevant_enc==1]$note_period_index)),
    "percentage of note periods validated by an procedure" = round((
       length(unique(output_procedures_exp[relevant_enc==1]$note_period_index))/
       length(unique(output[active_note==1]$note_period_index)))*100, digits=0),
    "number of notes validated by an procedure" = 
       length(unique(output_procedures_exp[relevant_enc==1]$note_id)),     
    "percentage of notes validated by an procedure" = round((
       length(unique(output_procedures_exp[relevant_enc==1]$note_id))/
       length(unique(output[active_note==1]$note_id)))*100, digits=0),

    "number of note periods validated by a diagnosis" = 
      length(unique(output_diagnoses_exp[relevant_enc==1]$note_period_index)),
    "percentage of note periods validated by a diagnosis" = round((
       length(unique(output_diagnoses_exp[relevant_enc==1]$note_period_index))/
       length(unique(output[active_note==1]$note_period_index)))*100, digits=0),
    "number of notes validated by a diagnosis" = 
       length(unique(output_diagnoses_exp[relevant_enc==1]$note_id)),     
    "percentage of notes validated by a diagnosis" = round((
       length(unique(output_diagnoses_exp[relevant_enc==1]$note_id))/
       length(unique(output[active_note==1]$note_id)))*100, digits=0),

    "number of note periods validated by a diagnosis and/or procedure and/or encounter" = 
      length(unique(c(output_diagnoses_exp[relevant_enc==1]$note_period_index,
        output_procedures_exp[relevant_enc==1]$note_period_index,
        output_encounters_exp[relevant_enc==1]$note_period_index))),
    "percentage of note periods validated by a diagnosis and/or procedure and/or encounter" = round((
       length(unique(c(output_diagnoses_exp[relevant_enc==1]$note_period_index,
        output_procedures_exp[relevant_enc==1]$note_period_index,
        output_encounters_exp[relevant_enc==1]$note_period_index)))/
       length(unique(output[active_note==1]$note_period_index)))*100, digits=0),
    "number of notes validated by a diagnosis and/or procedure and/or encounter" = 
       length(unique(c(output_diagnoses_exp[relevant_enc==1]$note_id,
        output_procedures_exp[relevant_enc==1]$note_id,
        output_encounters_exp[relevant_enc==1]$note_id))),     
    "percentage of notes validated by a diagnosis and/or procedure and/or encounter" = round((
       length(unique(c(output_diagnoses_exp[relevant_enc==1]$note_id,
        output_procedures_exp[relevant_enc==1]$note_id,
        output_encounters_exp[relevant_enc==1]$note_id)))/
       length(unique(output[active_note==1]$note_id)))*100, digits=0),
    "number of patients with a note that is validated by a diagnosis and/or 
    procedure and/or encounter" =
      length(unique(c(output_diagnoses_exp[relevant_enc==1]$empi,
        output_procedures_exp[relevant_enc==1]$empi,
        output_encounters_exp[relevant_enc==1]$empi))),  
    
    "number of ed encounters falling within the periods" =
      length(unique(output_ed_encounters_exp[relevant_enc==1]$encounter_number)),
    "number of patients with an ed encounter falling within the periods" =
      length(unique(output_ed_encounters_exp[relevant_enc==1]$empi)),
    
  
    "number of ed encounters falling within the periods by patients with a 
    diagnosis/procedure/encoutner at any time" =
      length(unique(output_ed_encounters_exp[relevant_enc==1 & empi %in%
      c(output_procedures_exp[!is.na(encounter_number) & 
      active_note==1]$empi,output_diagnoses_exp[!is.na(encounter_number)
      & active_note==1]$empi, output_encounters_exp[!is.na(encounter_number) & 
      active_note==1]$empi)]$encounter_number)),
    "number of patients with an ed encounters falling within the periods who 
    had a diagnosis/procedure/encoutner at any time" =
      length(unique(output_ed_encounters_exp[relevant_enc==1 & empi %in%
        c(output_procedures_exp[!is.na(encounter_number) & 
      active_note==1]$empi,output_diagnoses_exp[!is.na(encounter_number)
      & active_note==1]$empi, output_encounters_exp[!is.na(encounter_number) & 
      active_note==1]$empi)]$empi)),

    "number of ed encounters falling within the periods validated by a diagnosis 
     and/or procedure and/or encounter" =
      length(unique(output_ed_encounters_exp[relevant_enc==1 & note_period_index %in%
        c(output_procedures_exp[relevant_enc==1]$note_period_index, 
          output_diagnoses_exp[relevant_enc==1]$note_period_index,
          output_encounters_exp[relevant_enc==1]$note_period_index)]$encounter_number)),
    "number of patients with an ed encounters falling within the periods validated 
     by a diagnosis and/or procedure and/or encounter" =
      length(unique(output_ed_encounters_exp[relevant_enc==1 & note_period_index %in%
        c(output_procedures_exp[relevant_enc==1]$note_period_index, 
          output_diagnoses_exp[relevant_enc==1]$note_period_index,
          output_encounters_exp[relevant_enc==1]$note_period_index)]$empi))
  ) 

  targex_output_summary_stats_validation<-as.data.frame(t(as.data.frame(
    targex_output_summary_stats_validation)))
  targex_output_summary_stats_validation$V0 <- 
    c(rownames(targex_output_summary_stats_validation))
  colnames(targex_output_summary_stats_validation) <- 
    c("var", "var_name")
  targex_output_summary_stats_validation <- 
    as.data.table(targex_output_summary_stats_validation)

# SAVE - FINAL OUTPUT 1b - Summary stats at the note level
write.csv(targex_output_summary_stats_validation, file=paste0(output_data_folder,
  "targex_output_summary_stats_validation", filename,".csv"), row.names=F)
}

################################################################################
################ (F) OUTPUTS  II  - ACROSS PATINETS/TIME (AGGREGATE) ###########
################################################################################
### number of active notes per patient
plot1 <- ggplot(data=output[active_note==1], aes(x=active_notes_per_patient)) +
geom_histogram(binwidth=1, position="dodge")+
geom_vline(aes(xintercept=median(active_notes_per_patient)),colour="grey70")+
geom_text(aes(mean(active_notes_per_patient)+10,2000,label = 
paste("median number of \n notes/patient:\n", round(median(active_notes_per_patient),
digits=0))), size=1.5, colour = "grey50", family="URWHelvetica") + 
ggtitle("distribution of extracted notes across patients\n") +
ylab("number of patients\n") +
xlab("number of extracted notes/patient") +
theme(axis.text.x = element_text(size=rel(0.5), colour = "grey50", 
  family="URWHelvetica"),
  axis.text.y = element_text(size=rel(0.5), colour = "grey50", 
  family="URWHelvetica"), 
  axis.title.x = element_text(size=rel(0.5), colour = "grey50", 
  family="URWHelvetica"),
  axis.title.y = element_text(size=rel(0.5), colour = "grey50", 
  family="URWHelvetica"),
  plot.title = element_text(size = 8,colour = "black",family="Helvetica" ), 
  plot.margin= 
  unit(c(0.5, 0.5, 0, 0.5), "cm"))

################################################################################
### active notes across time (month windows) 
plot2 <- ggplot(data=output[active_note==1], aes(x=reorder(paste(lno_month, 
   lno_year),as.numeric(strptime(lno_date, "%Y-%m-%d"))),fill=factor(lno_month)))+
geom_bar()+
ylab("number of extracted notes\n") +
xlab("date of extracted notes") +
scale_x_discrete(labels=c(rep(c(""), length=5),"2010","", rep(c(""), length=5), 
  rep(c(""), length=5),"2011","", rep(c(""), length=5),
  rep(c(""), length=5),"2012","", rep(c(""), length=5),
  rep(c(""), length=5),"2013","", rep(c(""), length=5)))+
scale_fill_manual(values = rep(c(rgb(216,229,190,max=255),rgb(174,210,192,
  max=255),rgb(134,191,193,max=255),rgb(100,174,197,max=255),rgb(76,152,192,
  max=255),rgb(63,126,173,max=255),rgb(50,100,155,max=255),rgb(39,77,138,
  max=255),rgb(27,55,121,max=255),rgb(18,37,106,max=255),rgb(14,28,80,max=255),
  rgb(12,20,57,max=255)),3), 
labels=c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug","sep","oct",
  "nov","dec"),name="note month") +
ggtitle("distribution of extracted notes across time (months)\n") +
theme(axis.text.x = element_text(size=rel(0.5), colour = "grey50",
  family="URWHelvetica"),
  axis.text.y = element_text(size=rel(0.5), colour = "grey50",
  family="URWHelvetica"), 
  axis.title.x = element_text(size=rel(0.5), colour = "grey50",
  family="URWHelvetica"),
  axis.title.y = element_text(size=rel(0.5), colour = "grey50",
  family="URWHelvetica"),
  plot.title = element_text(size = 8,colour = "black",
  family="URWHelvetica"),
  legend.title = element_text(size = 5,colour = "grey50", 
  family="URWHelvetica"),
  legend.text = element_text(size = rel(0.5),colour = "grey50", 
  family="URWHelvetica"),  
  legend.position="bottom", legend.key.size=unit(0.2,"cm"), 
  plot.margin= unit(c(0.3, 0.5, 0.5, 0.5), "cm"))

# SAVE - FINAL OUTPUT 2 - Notes across time and patients (aggregate)
pdf(paste0(output_data_folder,"targex_output_distribution", 
  filename,".pdf"))
grid.arrange(plot1, plot2, ncol=1)
dev.off()

################################################################################
######################### (G) OUTPUTS  III  - CONCEPTS #########################
################################################################################
################################################################################
### concept drivers
concept_driver_table <- output[active_note==1, 
  c("empi", "note_id", "concept_count_mod", "lno_date", paste0(concept_list_expanded, 
  "_dummy_sum")), with=F]
concept_driver_table <- as.data.table(gather(concept_driver_table, 
  concept_name, concept_dummy, -empi, -note_id, -concept_count_mod, -lno_date))

concept_driver_table <- concept_driver_table[, .(note_count=length(
  unique(note_id[concept_dummy!=0])),note_count_unique=length(
  unique(note_id[concept_dummy!=0 & concept_count_mod==1]))),   
  by=c("concept_name")][, .(concept_name, note_count, note_perc=(note_count/
  length(unique(output[active_note!=0]$note_id)))*100, note_perc_unique=
  (note_count_unique/length(unique(output[active_note!=0]$note_id)))*100)]

concept_driver_table[, concept_name:=tolower(gsub("_dummy_sum", "", 
  concept_name))]
concept_driver_table <- as.data.table(rbind(
  setnames(concept_driver_table[, c(1,2,4), with=F], names(concept_driver_table[, 
  c(1,2,3), with=F])), concept_driver_table[, c(1,2,3), with=F]))

concept_driver_table[, cat:=c(rep("note_unique_perc", length(concept_list_expanded)), 
  rep("note_perc", length(concept_list_expanded)))]
concept_driver_table[, cat_level:=c(rep(1, length(concept_list_expanded)), 
  rep(2, length(concept_list_expanded)))]

colour <- ifelse(concept_list_expanded %in% concept_list, "black", "grey50")
plot1 <- ggplot(data=concept_driver_table[order(-note_perc)], aes(
  x=factor(concept_name), y=note_perc)) +
geom_bar(stat="identity", aes(fill=reorder(factor(cat), cat_level)), 
  position="dodge")+ 
coord_flip() +
ggtitle("concepts driving the note identification\n") +
ylab("\n% of notes incl. this concept/incl. only this concept") +
xlab("concept\n") +
scale_fill_manual(values=c(rgb(181,199,222,max=255), "steelblue")) +
theme(
  axis.text.x = element_text(size=rel(0.2), colour = "grey50",
  family="URWHelvetica"),
  axis.text.y = element_text(size=rel(0.2), colour = colour,
  family="URWHelvetica"), 
  axis.title.x = element_text(size=rel(0.5), colour = "grey50",
  family="URWHelvetica"),
  axis.title.y = element_text(size=rel(0.5), colour = "grey50",
  family="URWHelvetica"),
  plot.title = element_text(size = 8,colour = "black",family="URWHelvetica"),
  legend.position="none", axis.ticks = element_line(colour = 'grey50', 
  size = 0.1), axis.ticks.length = unit(0.01, "cm"))

################################################################################
### correlation between concept members
concept_correlation <- melt(cor(output[active_note==1,
  .SD,.SDcols=names(output)[names(output) %like% "[A-Z]" & names(output) 
  %like% "\\."]]))
concept_correlation <- as.data.table(concept_correlation)[order(-X1)]

concept_correlation[, c("X1", "X2"):=.(tolower(gsub("_dummy_sum", "", 
  X1)), tolower(gsub("_dummy_sum", "", X2)))]

plot2 <- ggplot(data=concept_correlation, aes(
  x=X1, y=X2, fill=value)) +
geom_tile(colour="white")+ 
ggtitle("correlation between concepts\n") +
 scale_fill_gradient(low = "white",high = "steelblue") + 
theme_grey(base_size=10) + labs(x = "",y = "") + 
scale_x_discrete(expand= c(0,0)) +
scale_y_discrete(expand = c(0,0)) + 
theme(
  axis.text.x = element_text(size=rel(0.4), colour = "grey50",
  family="URWHelvetica", angle=90, hjust=1),
  axis.text.y = element_text(size=rel(0.4), colour = "grey50",
  family="URWHelvetica"), 
  axis.title.x = element_text(size=rel(0.5), colour = "grey50",
  family="URWHelvetica"),
  axis.title.y = element_text(size=rel(0.5), colour = "grey50",
  family="URWHelvetica"),
  plot.title = element_text(size = 8,colour = "black",family="URWHelvetica"),
  legend.position="none", axis.ticks = element_line(colour = 'grey50', 
  size = 0.1), axis.ticks.length = unit(0.01, "cm"), 
  plot.margin= unit(c(0.3, 3, 0.5, 2), "cm"))


# SAVE - FINAL OUTPUT 3 - Notes across concepts
pdf(paste0(output_data_folder,"targex_output_concepts", 
  filename,".pdf"))
grid.arrange(plot1, plot2, ncol=1, heights=c(15,20),widths=10)
dev.off()

################################################################################
###################### (H) OUTPUTS  IV  - PATIENT LEVEL ########################
################################################################################
### notes across time for subset of patients & illustrate subseeting processs if 
### implemented

# subset relevant data & patients
sample <- empi_index_expanded
sample_size <- length(sample$empi)

output_patient_time <- output[, .(empi, lno_year_week, active_note,
  active_note_broad)]

# generate weekly dummies
for (var in paste0(rep(c("2010-","2011-","2012-","2013-"),each=52),
  rep(sprintf("%02d", 1:52),4),"-1")) {
  output_patient_time[lno_year_week==var, paste0("D_",var) :=.(1)]
  output_patient_time[lno_year_week!=var, paste0("D_",var):=.(0)]
  output_patient_time[lno_year_week==var & active_note==1, 
    paste0("DA_",var) :=.(1)]
  output_patient_time[lno_year_week!=var & active_note==1,
    paste0("DA_",var):=.(0)]
  output_patient_time[active_note==0, paste0("DA_",var):=.(0)]
  output_patient_time[lno_year_week==var & active_note_broad==1, 
    paste0("DAB_",var) :=.(1)]
  output_patient_time[lno_year_week!=var & active_note_broad==1, 
    paste0("DAB_",var):=.(0)]
  output_patient_time[active_note_broad==0, paste0("DAB_",var):=.(0)]
}

output_patient_time[, c(paste0("D_",rep(c("2010-","2011-","2012-","2013-"),
  each=52),rep(sprintf("%02d", 1:52),4),"-1")) := lapply(.SD, function(x) 
  sum(as.numeric(x))), by=c("empi"), .SDcols=paste0("D_",rep(c("2010-","2011-",
  "2012-","2013-"),each=52),rep(sprintf("%02d", 1:52),4), "-1")]

output_patient_time[, c(paste0("DA_",rep(c("2010-","2011-","2012-",
  "2013-"),each=52),rep(sprintf("%02d", 1:52),4), "-1")) := lapply(.SD, 
  function(x) sum(as.numeric(x))), by=c("empi"), .SDcols=paste0("DA_",
  rep(c("2010-","2011-","2012-","2013-"),each=52),rep(sprintf("%02d", 1:52),4), 
  "-1")]

output_patient_time[, c(paste0("DAB_",rep(c("2010-","2011-","2012-",
  "2013-"),each=52),rep(sprintf("%02d", 1:52),4), "-1")) := lapply(.SD, 
  function(x) sum(as.numeric(x))), by=c("empi"), .SDcols=paste0("DAB_",
  rep(c("2010-","2011-","2012-","2013-"),each=52),rep(sprintf("%02d", 1:52),4), 
  "-1")]

output_patient_time <- unique(output_patient_time, by=c("empi"))

output_patient_time <- as.data.table(gather(output_patient_time,year_week,
  year_week_dummy, -lno_year_week, -empi, - active_note, - active_note_broad))

output_patient_time[year_week %like% "D_", category:="all"]
output_patient_time[year_week %like% "DA_", category:="active"]
output_patient_time[year_week %like% "DAB_", category:="active_broad"]

output_patient_time[, year_week:=as.IDate(gsub("D_|DA_|DAB_", "", year_week),
  "%Y-%W-%u")]

output_patient_time[year_week_dummy==0, year_week_dummy:=NA]

# subset to relevant patients
output_patient_time_sample <- output_patient_time[empi %in% sample$empi]

output_patient_time_sample <- output_patient_time_sample[order(empi, year_week)]
output_patient_time_sample[category=="all", patient_index:=rep(seq(1,1+
  (sample_size*0.05)-0.05, by=0.05), each=208)]
output_patient_time_sample[category=="active", patient_index:=rep(seq(1,1+
  (sample_size*0.05)-0.05, by=0.05), each=208)]
output_patient_time_sample[category=="active_broad", patient_index:=
  rep(seq(1,1+(sample_size*0.05)-0.05, by=0.05), each=208)]

# XXX TO-DO: Remove need for this gap year fill
fill_data_leap_year_2012<- output_patient_time_sample[ 
  year_week=="2012-10-01"]
fill_data_leap_year_2012[,
  c("year_week", "year_week_dummy"):=.(as.IDate("2012-53-1", "%Y-%W-%u"), NA)]
output_patient_time_sample <- as.data.table(rbind(output_patient_time_sample, 
  fill_data_leap_year_2012))


# SAVE - DATA OUTPUT 6 - Notes across time at the patient level - 
# subset of patients
write.csv(output_patient_time, file=paste0(modified_data_folder, 
  "targex_output_patient_time", filename, ".csv"), row.names=F)
write.csv(output_patient_time_sample, file=paste0(modified_data_folder, 
  "targex_output_patient_time_sample", filename, ".csv"), row.names=F)

# plot main
plot1 <- ggplot(data=output_patient_time_sample[category=="all"], 
  aes(x=as.numeric(strptime(year_week,"%Y-%m-%d")), y=patient_index))+
  geom_tile(aes(height=0.05, fill=year_week_dummy), colour="white") +
xlab("week")  +
ylab(paste0("patients (subset of ",sample_size, " patients)")) +
ggtitle("extracted notes across time for a representative patient subset (all)\n") + 
scale_fill_gradient2(low="orange", high="red", na.value=rgb(245,245,245, 
    alpha=100, max=255))  +
scale_x_continuous(breaks = seq(min(as.numeric(
  strptime(output_patient_time_sample[category=="all"]$year_week,"%Y-%m-%d"))), 
  max(as.numeric(strptime(output_patient_time_sample[category=="all"]$year_week,
  "%Y-%m-%d"))), length.out=208), limits=c(min(as.numeric(
  strptime(output_patient_time_sample[category=="all"]$year_week,"%Y-%m-%d"))), 
  max(as.numeric(strptime(output_patient_time_sample[category=="all"]$year_week,
  "%Y-%m-%d")))), labels=c(rep("",26), "2010",rep("",25), rep("",26),"2011", 
  rep("",25), rep("",26), "2012",rep("",25), rep("",26), "2013", rep("",25)), 
  expand=c(0,0))+
 geom_segment(x=c(as.numeric(strptime("2011-01-01","%Y-%m-%d"))), 
    xend=c(as.numeric(strptime("2011-01-01","%Y-%m-%d"))),
    y=1, yend=1+(sample_size*0.05)-0.05, colour="grey80", size=0.) +
 geom_segment(x=c(as.numeric(strptime("2012-01-01","%Y-%m-%d"))), 
   xend=c(as.numeric(strptime("2012-01-01","%Y-%m-%d"))),
   y=1, yend=1+(sample_size*0.05)-0.05, colour="grey80", size=0.) +
  geom_segment(x=c(as.numeric(strptime("2013-01-01","%Y-%m-%d"))), 
   xend=c(as.numeric(strptime("2013-01-01","%Y-%m-%d"))),
   y=1, yend=1+(sample_size*0.05)-0.05, colour="grey80", size=0.) +
  scale_y_continuous(expand=c(0,0))+
theme(  
  axis.title.y=element_text(colour="grey50", size=rel(0.4)),
  axis.title.x=element_text(colour="grey50", size=rel(0.4)),
  axis.text.x = element_text(size=rel(0.3), colour = "grey50"), 
  axis.text.y = element_blank(), 
  plot.title=element_text(colour="grey50", size=6),
  panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
  panel.background = element_rect(fill ="white"),
   plot.background = element_rect(fill ="white"),
  legend.title=element_text(colour="grey50", size=rel(0.8), face="plain"),
  legend.text=element_text(colour="grey50", size=rel(0.7)), 
  axis.ticks = element_blank(), legend.position="none",plot.margin= 
  unit(c(0.5, 0.5, 0, 0.5), "cm"))


# additional plots if subsetting is implemented
if (!is.na(time_frame_frequency_weeks)) {


  output_data <- output[output_patient_time_sample[category=="all"], on=c("empi")]

  plot2 <- ggplot(data=output_patient_time_sample[category=="all"], 
    aes(x=as.numeric(strptime(year_week,"%Y-%m-%d")), y=patient_index))+
  geom_tile(aes(height=0.05, fill=year_week_dummy), colour="white") +
  scale_fill_gradient2(low="orange", high="red", na.value=rgb(245,245,245, 
    alpha=100, max=255))  +
 xlab("week")  +
 ylab(paste0("patients (subset of ",sample_size, " patients)")) +
 ggtitle(
   "extracted notes across time for a representative patient subset (all and subset)\n") + 
 geom_rect(data=output_patient_time_sample[category=="active_broad"][!
    year_week_dummy==0], aes(xmin=as.numeric(strptime(year_week,"%Y-%m-%d")) - 
    303861, xmax=as.numeric(strptime(year_week,"%Y-%m-%d")) +303861, 
    ymin=patient_index-0.025 , ymax=patient_index+0.025 ), colour="black", fill=NA,
    size=0.1) +
 geom_rect(data=output_data[active_note_broad==1 & note_period %in% c(0,3)], 
  	aes(xmin=as.numeric(strptime(lno_year_week,"%Y-%W-%u")) - 
    303861, xmax=as.numeric(strptime(lno_year_week,"%Y-%W-%u")) +303861, 
    ymin=patient_index-0.025 , ymax=patient_index+0.025 ), colour="black", 
    fill="black", size=0.1) +
 geom_point(data=output_patient_time_sample[category=="active"][!
    year_week_dummy==0], aes(x=as.numeric(strptime(year_week,"%Y-%m-%d")), 
    ymin=patient_index), fill="darkred",colour="darkred", shape=5,size=0.1) +
 scale_x_continuous(breaks = seq(min(as.numeric(
    strptime(output_patient_time_sample[category=="all"]$year_week,"%Y-%m-%d"))), 
    max(as.numeric(strptime(output_patient_time_sample[category=="all"]$year_week,
    "%Y-%m-%d"))), length.out=208), limits=c(min(as.numeric(
    strptime(output_patient_time_sample[category=="all"]$year_week,"%Y-%m-%d"))), 
    max(as.numeric(strptime(output_patient_time_sample[category=="all"]$year_week,
    "%Y-%m-%d")))), labels=c(rep("",26), "2010",rep("",25), rep("",26),"2011", 
    rep("",25), rep("",26), "2012",rep("",25), rep("",26), "2013", rep("",25)), 
    expand=c(0,0))+
 scale_y_continuous(expand=c(0,0))+
 geom_segment(x=c(as.numeric(strptime("2011-01-01","%Y-%m-%d"))), 
    xend=c(as.numeric(strptime("2011-01-01","%Y-%m-%d"))),
    y=1, yend=1+(sample_size*0.05)-0.05, colour="grey80", size=0.1) +
 geom_segment(x=c(as.numeric(strptime("2012-01-01","%Y-%m-%d"))), 
   xend=c(as.numeric(strptime("2012-01-01","%Y-%m-%d"))),
   y=1, yend=1+(sample_size*0.05)-0.05, colour="grey80", size=0.1) +
  geom_segment(x=c(as.numeric(strptime("2013-01-01","%Y-%m-%d"))), 
   xend=c(as.numeric(strptime("2013-01-01","%Y-%m-%d"))),
   y=1, yend=1+(sample_size*0.05)-0.05, colour="grey80", size=0.1) +
  theme(  
   axis.title.y=element_text(colour="grey50", size=rel(0.4)),
  axis.title.x=element_text(colour="grey50", size=rel(0.4)),
  axis.text.x = element_text(size=rel(0.3), colour = "grey50"), 
  axis.text.y = element_blank(), 
  plot.title=element_text(colour="grey50", size=6),
  panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
  panel.background = element_rect(fill ="white"),
   plot.background = element_rect(fill ="white"),
  legend.title=element_text(colour="grey50", size=rel(0.8), face="plain"),
  legend.text=element_text(colour="grey50", size=rel(0.7)), 
  axis.ticks = element_blank(), legend.position="none",plot.margin= 
  unit(c(0.5, 0.5, 0, 0.5), "cm"))


  plot3 <- ggplot(data=output_patient_time_sample[category=="active"], 
  aes(x=as.numeric(strptime(year_week,"%Y-%m-%d")), y=patient_index))+
  geom_tile(aes(height=0.05, fill=year_week_dummy), colour="white") +
xlab("week")  +
ylab(paste0("patients (subset of ",sample_size, " patients)")) +
ggtitle(
  "extracted notes across time for a representative patient subset (subset)\n") + 
scale_fill_gradient2(low="orange", high="red", na.value=rgb(245,245,245, 
    alpha=100, max=255))  +
scale_x_continuous(breaks = seq(min(as.numeric(
  strptime(output_patient_time_sample[category=="all"]$year_week,"%Y-%m-%d"))), 
  max(as.numeric(strptime(output_patient_time_sample[category=="all"]$year_week,
  "%Y-%m-%d"))), length.out=208), limits=c(min(as.numeric(
  strptime(output_patient_time_sample[category=="all"]$year_week,"%Y-%m-%d"))), 
  max(as.numeric(strptime(output_patient_time_sample[category=="all"]$year_week,
  "%Y-%m-%d")))), labels=c(rep("",26), "2010",rep("",25), rep("",26),"2011", 
  rep("",25), rep("",26), "2012",rep("",25), rep("",26), "2013", rep("",25)), 
  expand=c(0,0))+
 geom_segment(x=c(as.numeric(strptime("2011-01-01","%Y-%m-%d"))), 
    xend=c(as.numeric(strptime("2011-01-01","%Y-%m-%d"))),
    y=1, yend=1+(sample_size*0.05)-0.05, colour="grey80", size=0.) +
 geom_segment(x=c(as.numeric(strptime("2012-01-01","%Y-%m-%d"))), 
   xend=c(as.numeric(strptime("2012-01-01","%Y-%m-%d"))),
   y=1, yend=1+(sample_size*0.05)-0.05, colour="grey80", size=0.) +
  geom_segment(x=c(as.numeric(strptime("2013-01-01","%Y-%m-%d"))), 
   xend=c(as.numeric(strptime("2013-01-01","%Y-%m-%d"))),
   y=1, yend=1+(sample_size*0.05)-0.05, colour="grey80", size=0.) +
  scale_y_continuous(expand=c(0,0))+
theme(  
  axis.title.y=element_text(colour="grey50", size=rel(0.4)),
  axis.title.x=element_text(colour="grey50", size=rel(0.4)),
  axis.text.x = element_text(size=rel(0.3), colour = "grey50"), 
  axis.text.y = element_blank(), 
  plot.title=element_text(colour="grey50", size=6),
  panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
  panel.background = element_rect(fill ="white"),
   plot.background = element_rect(fill ="white"),
  legend.title=element_text(colour="grey50", size=rel(0.8), face="plain"),
  legend.text=element_text(colour="grey50", size=rel(0.7)), 
  axis.ticks = element_blank(), legend.position="none",plot.margin= 
  unit(c(0.5, 0.5, 0, 0.5), "cm"))

  # SAVE FINAL OUTPUT 4 - Notes across time for patients
  pdf(paste0(output_data_folder, "targex_output_patient_time", filename, 
  	".pdf"))
  grid.arrange(plot1,plot2,plot3, heights=c(2,2,2), widths=20, ncol=1)
  dev.off()

} else {

  # SAVE - FINAL OUTPUT 4 - Notes across time for patients
  ggsave(filename=paste0(output_data_folder, "targex_output_patient_time", 
  filename, ".pdf"),plot=plot1, height = 3, width = 10)

}

################################################################################
### notes across time for patients - concepts, notes and ed encoutners vs.
### concepts, notes and encoutners, diagnsoes, procedures - mark relevant 
### encoutners if subsetting is implemented

# subset relevant data & patients
sample <- empi_index_high_frequency
sample_size <- length(unique(sample$empi))
variable_number <- ifelse(variables_to_graph=="concept_list_expanded", 
  length(concept_list_expanded), length(concept_list))

if (variables_to_graph=="concept_list_expanded") {
output_patient_time_concept <- output[, .SD, .SDcols=union(c("empi", 
   "lno_year_week"),names(output)[names(output) %like% "[A-Z]" & 
   names(output) %like% "_dummy_sum"])][empi %in% sample$empi]
} else {
output_patient_time_concept <- output[, .SD, .SDcols=union(c("empi", 
  "lno_year_week"),paste0(concept_list,"_dummy_sum"))][
   empi %in% sample$empi]
}

# generate weekly dummies 
for (var in paste0(rep(c("2010-","2011-","2012-","2013-"),each=52),
  rep(sprintf("%02d", 1:52),4),"-1")) {
  output_patient_time_concept[lno_year_week==var, paste0("D_",var) :=.(1)]
  output_patient_time_concept[lno_year_week!=var, paste0("D_",var):=.(0)]
}

for (var in c(names(output_patient_time_concept)[names(output_patient_time_concept) 
  %like% "[A-Z]" & !names(output_patient_time_concept) %like% "---|D_" ])) {
  for (date_var in c(paste0("D_",rep(c("2010-","2011-","2012-","2013-"),
  each=52),rep(sprintf("%02d", 1:52),4), "-1") )) {
    output_patient_time_concept[,c(paste0(var, "---", date_var)):=do.call(
    paste,c(.SD, sep="---")), .SDcols=c(var, date_var)]
  }
}

output_patient_time_concept[, 
c(names(output_patient_time_concept)[names(output_patient_time_concept) 
  %like% "[A-Z]" & !names(output_patient_time_concept) %like% "---|D_" ], 
  paste0("D_",rep(c("2010-","2011-","2012-","2013-"),each=52),
  rep(sprintf("%02d", 1:52),4),"-1")):=NULL]

output_patient_time_concept <- as.data.table(gather(output_patient_time_concept, 
  name, name_dummy, -empi, -lno_year_week))

output_patient_time_concept[,c("concept_name", "year_week"):=
  tstrsplit(name, "---")][, c("name"):=NULL]
output_patient_time_concept[,c("concept_name_dummy", "year_week_dummy"):=
  tstrsplit(name_dummy, "---")][,c("name_dummy"):=NULL]

output_patient_time_concept[year_week_dummy==0, concept_name_dummy:=0]

output_patient_time_concept[,paste0("year_week_dummy","_new"):=
  sum(as.numeric(year_week_dummy)), by=c("empi", "year_week", 
  "concept_name")]
output_patient_time_concept[,paste0("concept_name_dummy","_new"):=
  sum(as.numeric(concept_name_dummy)), by=c("empi", "year_week", 
  "concept_name")]
output_patient_time_concept[,c("concept_name_dummy",
  "year_week_dummy"):=NULL]

output_patient_time_concept <- unique(output_patient_time_concept, 
  by=c("empi", "year_week", "concept_name"))

output_patient_time_concept[, year_week:=as.IDate(gsub("D_", 
  "", year_week),"%Y-%W-%u")]

output_patient_time_concept <- output_patient_time_concept[
  order(empi, -concept_name, year_week)]

  # XXX TO-DO: Think of a better way of ensuring that concepts and 
  # cocnept members are appropriately sorted

  if (variables_to_graph=="concept_list_expanded") {
    temp_concept_list <- as.data.table(gsub("_dummy_sum", "", concept_list))[, 
      concept_id:=1:.N]
    temp_concept_list_expanded <- as.data.table(unique(
      output_patient_time_concept$concept_name))[, output_id:=1:.N]

    for (i in 1:nrow(temp_concept_list)) {
      temp_concept_list_expanded[gsub("\\..+|_dummy_sum", "", V1) %in% 
        temp_concept_list[i]$V1 , concept_exp_id:=temp_concept_list[i]$concept_id]
    }

    temp_concept_list_expanded <- temp_concept_list_expanded[order(concept_exp_id, 
      V1)][, concept_exp_id_sub:=1:.N, by=c("concept_exp_id")]
    
    temp_concept_list_expanded <- temp_concept_list_expanded[order(concept_exp_id, 
      concept_exp_id_sub)][, concept_exp_id_comb:=1:.N][order(output_id)]

    output_patient_time_concept  <- output_patient_time_concept[, concept_name_id:=
      temp_concept_list_expanded$concept_exp_id_comb, by=c("empi", "year_week")]

    output_patient_time_concept <- output_patient_time_concept[order(empi, 
      concept_name_id, year_week)]
  } 


output_patient_time_concept[,id:=rep(seq(1,(1+
  (((sample_size*variable_number)-1)*0.05)),by=0.05),each=208)]

# XXX TO-DO: Remove need for this gap year fill
fill_data_leap_year_2012<- output_patient_time_concept[
  year_week=="2012-10-01" ]
fill_data_leap_year_2012[,
  c("year_week", "concept_name_dummy_new" ):=.(as.IDate("2012-53-1", 
  "%Y-%W-%u"), 0)]
output_patient_time_concept <- as.data.table(rbind(output_patient_time_concept, 
  fill_data_leap_year_2012))

# SAVE - DATA OUTPUT 7 - Notes across time/concepts at the patient level - 
# subset of patients
write.csv(output_patient_time_concept, file=paste0(modified_data_folder, 
  "targex_output_patient_time_concept", filename, ".csv", row.names=F))

# plot main
for (i in seq(1,(1+(((sample_size*variable_number)-1)*0.05)),
  by=variable_number*0.05)) {
  
  data_plot<- output_patient_time_concept[id %in% seq(i,i+(0.05*variable_number)
    -0.05, by=0.05)]

  data_plot[, concept_name_dummy_new_new:=concept_name_dummy_new][
    year_week_dummy_new==0, concept_name_dummy_new_new:=0]

  data_plot[concept_name_dummy_new_new==0, concept_name_dummy_new_new:=NA]

  notes_plot_all <- output_patient_time[category=="all"][empi 
    %in% data_plot$empi]
  notes_plot_all[, id:=i+((variable_number+0.5)*0.05)]
  
  # XXX TO-DO: Remove need for this gap year fill
  fill_data_leap_year_2012<- notes_plot_all[
  year_week=="2012-10-01" ]
  fill_data_leap_year_2012[,c("year_week", "year_week_dummy" ):=
  .(as.IDate("2012-53-1", "%Y-%W-%u"), 0)]
  notes_plot_all <- as.data.table(rbind(notes_plot_all, 
  fill_data_leap_year_2012))

  notes_plot_all[year_week_dummy==0, year_week_dummy:=NA]

  notes_plot_active<- output_patient_time[category=="active"][
    empi %in% data_plot$empi][!year_week_dummy==0]
  notes_plot_active[, id:=i+((variable_number+0.5)*0.05)]

  # XXX TO-DO: Remove need for this gap year fill
  fill_data_leap_year_2012<- notes_plot_active[
  year_week=="2012-10-01" ]
  fill_data_leap_year_2012[,c("year_week", "year_week_dummy" ):=
  .(as.IDate("2012-53-1", "%Y-%W-%u"), 0)]
  notes_plot_active <- as.data.table(rbind(notes_plot_active, 
  fill_data_leap_year_2012))

  notes_plot_active[year_week_dummy==0, year_week_dummy:=NA]

 all_weeks <- as.data.table(c(paste0(rep(c("2010-","2011-"),
  	each=52),rep(sprintf("%02d", 1:52),2),"-1"), 
  	paste0(rep(c("2012-"),each=53),rep(sprintf("%02d", 1:53),1),"-1"), 
    paste0(rep(c("2013-"),each=52),rep(sprintf("%02d", 1:52),1),"-1")))

  encounters_plot <- unique(output_encounters_coll[!(cat=="active" & 
  enc_count==0)][order(cat)], by=c("empi", "admit_year_week")) [empi %in% 
    data_plot$empi][all_weeks, on=c(admit_year_week="V1"), nomatch=NA][, 
    admit_year_week:=as.IDate(admit_year_week,"%Y-%W-%u")][is.na(enc_count), 
    enc_count:=0]

  encounters_plot[enc_count==0, enc_count:=NA]

  ed_encounters_plot <- unique(output_ed_encounters_coll[!(cat=="active" & 
  ed_enc_count==0)][order(cat)], by=c("empi", "ed_year_week")) [empi %in% 
    data_plot$empi][all_weeks, on=c(ed_year_week="V1"), nomatch=NA][, 
    ed_year_week:=as.IDate(ed_year_week,"%Y-%W-%u")][is.na(ed_enc_count), 
    ed_enc_count:=0]

  ed_encounters_plot[ed_enc_count==0, ed_enc_count:=NA]


  ed_encounters_plot_act <- unique(output_ed_encounters_coll[
    !(cat=="active" & ed_enc_count==0)][order(cat)], by=c("empi", "ed_year_week")) [empi %in% 
    data_plot$empi][all_weeks, on=c(ed_year_week="V1"), nomatch=NA][, 
    ed_year_week:=as.IDate(ed_year_week,"%Y-%W-%u")][is.na(ed_enc_count), 
    ed_enc_count:=0][cat=="active"]

  ed_encounters_plot_act[ed_enc_count==0, ed_enc_count:=NA]

  diagnoses_plot <- unique(output_diagnoses_coll[!(cat=="active" & 
  diagnosis_count==0)][order(cat)], by=c("empi", "diagnosis_year_week")) [empi %in% 
    data_plot$empi][all_weeks, on=c(diagnosis_year_week="V1"), nomatch=NA][, 
    diagnosis_year_week:=as.IDate(diagnosis_year_week,"%Y-%W-%u")][
    is.na(diagnosis_count), diagnosis_count:=0]

  diagnoses_plot[diagnosis_count==0, diagnosis_count:=NA]

  procedures_plot <- unique(output_procedures_coll[!(cat=="active" & 
  procedure_count==0)][order(cat)], by=c("empi", "procedure_year_week")) [empi %in% 
    data_plot$empi][all_weeks, on=c(procedure_year_week="V1"), nomatch=NA][, 
    procedure_year_week:=as.IDate(procedure_year_week,"%Y-%W-%u")][
    is.na(procedure_count), procedure_count:=0]
 
  procedures_plot[procedure_count==0, procedure_count:=NA]

  colours_1 <- c(ifelse(unique(data_plot$concept_name) %like% "\\.", 
  	"grey50", "black"), "black", "black")

  colours_2 <- c(ifelse(unique(data_plot$concept_name) %like% "\\.", 
  	"grey50", "black"), "black", "black", "black", "black")


  plot1 <- ggplot(data=data_plot, 
  aes(x=as.numeric(strptime(year_week,"%Y-%m-%d")), y=id))+
  geom_tile(aes(height=0.05, fill=concept_name_dummy_new_new), colour="white") +
  geom_tile(data=notes_plot_all, aes(height=0.05, fill=year_week_dummy), 
   colour="white") +
  geom_tile(data=ed_encounters_plot, aes(x=as.numeric(strptime(ed_year_week,
    "%Y-%m-%d")), y=i+((variable_number+0.5)*0.05)+0.075, 
     height=0.05, fill=ed_enc_count), colour="white") +
  scale_fill_gradient2(low="orange", high="red", na.value=rgb(245,245,245, 
    alpha=100, max=255))  + 
  scale_y_continuous(breaks=c(unique(data_plot$id), i+((variable_number+0.5)*0.05),
  	i+((variable_number+0.5)*0.05)+0.075), 
    labels=c(tolower(gsub("_dummy_sum", "", unique(data_plot$concept_name))), 
    "extracted notes", "ed encounters"),
    expand=c(0,0))+
  xlab(paste0("week (patient - empi: ", unique(data_plot$empi), ")"))  +
  ylab("")+
  ggtitle("extracted concepts/notes and ed encounters across time (patient level)\n") + 
  scale_x_continuous(breaks = seq(min(as.numeric(
    strptime(data_plot$year_week,"%Y-%m-%d"))), 
    max(as.numeric(strptime(data_plot$year_week,
    "%Y-%m-%d"))), length.out=208), limits=c(min(as.numeric(
    strptime(data_plot$year_week,"%Y-%m-%d"))), 
    max(as.numeric(strptime(data_plot$year_week,
    "%Y-%m-%d")))), labels=c(rep("",26), "2010",rep("",25), rep("",26),"2011", 
    rep("",25), rep("",26), "2012",rep("",25), rep("",26), "2013", rep("",25)),
    expand=c(0,0))+
 geom_segment(x=c(as.numeric(strptime("2011-01-01","%Y-%m-%d"))), 
    xend=c(as.numeric(strptime("2011-01-01","%Y-%m-%d"))),
    y=i, yend=i+(0.05*variable_number)-0.05, colour="grey80", size=0.1) +
 geom_segment(x=c(as.numeric(strptime("2012-01-01","%Y-%m-%d"))), 
   xend=c(as.numeric(strptime("2012-01-01","%Y-%m-%d"))),
   y=i, yend=i+(0.05*variable_number)-0.05, colour="grey80", size=0.1) +
  geom_segment(x=c(as.numeric(strptime("2013-01-01","%Y-%m-%d"))), 
   xend=c(as.numeric(strptime("2013-01-01","%Y-%m-%d"))),
   y=i, yend=i+(0.05*variable_number)-0.05, colour="grey80", size=0.1) +
  theme(  
  axis.title.y=element_text(colour="grey50", size=rel(0.4), 
  	family="URWHelvetica"),
  axis.title.x=element_text(colour="grey50", size=rel(0.4), 
  	family="URWHelvetica"),
  axis.text.x = element_text(size=rel(0.3), colour = "grey50", 
  	family="URWHelvetica", hjust=1), 
  axis.text.y = element_text(size=rel(0.3), colour =  colours_1,hjust=1), 
  plot.title=element_text(colour="grey50", size=6, family="URWHelvetica"),
  panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
  panel.background = element_rect(fill ="white"),
   plot.background = element_rect(fill ="white"),
  legend.title=element_text(colour="grey50", size=rel(0.8), face="plain"),
  legend.text=element_text(colour="grey50", size=rel(0.7)), 
  axis.ticks = element_blank(),
  legend.position="none",legend.key.size=unit(0.2,"cm"), plot.margin= 
  unit(c(0.5, 1.2, 0, 0.2), "cm")) 

  plot2 <- ggplot(data=data_plot, 
  aes(x=as.numeric(strptime(year_week,"%Y-%m-%d")), y=id))+
  geom_tile(aes(height=0.05, fill=concept_name_dummy_new_new), colour="white") +
  geom_tile(data=notes_plot_all, aes(height=0.05, fill=year_week_dummy), 
   colour="white") +
  geom_tile(data=encounters_plot, aes(x=as.numeric(strptime(admit_year_week,
    "%Y-%m-%d")), y=i+((variable_number+0.5)*0.05)+0.075, 
    height=0.05, fill=enc_count), colour="white") +
  geom_tile(data=diagnoses_plot, aes(x=as.numeric(strptime(diagnosis_year_week,
    "%Y-%m-%d")), y=i+((variable_number+0.5)*0.05)+0.15, 
    height=0.05, fill=diagnosis_count), colour="white") +
  geom_tile(data=procedures_plot, aes(x=as.numeric(strptime(procedure_year_week,
    "%Y-%m-%d")), y=i+((variable_number+0.5)*0.05)+0.225, 
    height=0.05, fill=procedure_count), colour="white") +
  scale_fill_gradient2(low="orange", high="red", na.value=rgb(245,245,245, 
    alpha=100, max=255))  + 
  scale_y_continuous(breaks=c(unique(data_plot$id), i+((variable_number+0.5)*0.05),
  	i+((variable_number+0.5)*0.05)+0.075, i+((variable_number+0.5)*0.05) + 0.15,
  	i+((variable_number+0.5)*0.05)+0.225), 
    labels=c(tolower(gsub("_dummy_sum", "", unique(data_plot$concept_name))), 
    "extracted notes", "rel. encounters", "rel. diagnoses", "rel. procedures" ),
    expand=c(0,0))+
  xlab(paste0("week (patient - empi: ", unique(data_plot$empi), ")"))  +
  ylab("")+
  ggtitle("extracted concepts/notes and ed encounters across time (patient level)\n") + 
  scale_x_continuous(breaks = seq(min(as.numeric(
    strptime(data_plot$year_week,"%Y-%m-%d"))), 
    max(as.numeric(strptime(data_plot$year_week,
    "%Y-%m-%d"))), length.out=208), limits=c(min(as.numeric(
    strptime(data_plot$year_week,"%Y-%m-%d"))), 
    max(as.numeric(strptime(data_plot$year_week,
    "%Y-%m-%d")))), labels=c(rep("",26), "2010",rep("",25), rep("",26),"2011", 
    rep("",25), rep("",26), "2012",rep("",25), rep("",26), "2013", rep("",25)),
    expand=c(0,0))+
 geom_segment(x=c(as.numeric(strptime("2011-01-01","%Y-%m-%d"))), 
    xend=c(as.numeric(strptime("2011-01-01","%Y-%m-%d"))),
    y=i, yend=i+(0.05*variable_number)-0.05, colour="grey80", size=0.1) +
 geom_segment(x=c(as.numeric(strptime("2012-01-01","%Y-%m-%d"))), 
   xend=c(as.numeric(strptime("2012-01-01","%Y-%m-%d"))),
   y=i, yend=i+(0.05*variable_number)-0.05, colour="grey80", size=0.1) +
  geom_segment(x=c(as.numeric(strptime("2013-01-01","%Y-%m-%d"))), 
   xend=c(as.numeric(strptime("2013-01-01","%Y-%m-%d"))),
   y=i, yend=i+(0.05*variable_number)-0.05, colour="grey80", size=0.1) +
  theme(  
  axis.title.y=element_text(colour="grey50", size=rel(0.4), 
  	family="URWHelvetica"),
  axis.title.x=element_text(colour="grey50", size=rel(0.4), 
  	family="URWHelvetica"),
  axis.text.x = element_text(size=rel(0.3), colour = "grey50", 
  	family="URWHelvetica", hjust=1), 
  axis.text.y = element_text(size=rel(0.3), colour =  colours_2,hjust=1), 
  plot.title=element_text(colour="grey50", size=6, family="URWHelvetica"),
  panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
  panel.background = element_rect(fill ="white"),
   plot.background = element_rect(fill ="white"),
  legend.title=element_text(colour="grey50", size=rel(0.8), face="plain"),
  legend.text=element_text(colour="grey50", size=rel(0.7)), 
  axis.ticks = element_blank(),
  legend.position="none",legend.key.size=unit(0.2,"cm"), plot.margin= 
  unit(c(0.5, 1.2, 0.5, 0.2), "cm"))  

  if (!is.na(time_frame_frequency_weeks)) {
  	plot1 <- plot1 +
    geom_rect(data=notes_plot_active[!
    year_week_dummy==0], aes(xmin=as.numeric(strptime(year_week,"%Y-%m-%d")) - 
    303861, xmax=as.numeric(strptime(year_week,"%Y-%m-%d")) +303861, 
    ymin=id-0.025 , ymax=id+0.025 ), colour="grey50", fill=NA,
    size=0.3) 

    if (nrow(ed_encounters_plot_act)!=0) {
      plot1 <- plot1 + geom_point(data=ed_encounters_plot_act, 
        aes(x=as.numeric(strptime(ed_year_week,
        "%Y-%m-%d")), y=i+((variable_number+0.5)*0.05)+0.075), colour="black", shape=5,
        size=0.1) 
    }

  	plot2 <- plot2 +
    geom_rect(data=notes_plot_active[!
    year_week_dummy==0], aes(xmin=as.numeric(strptime(year_week,"%Y-%m-%d")) - 
    303861, xmax=as.numeric(strptime(year_week,"%Y-%m-%d")) +303861, 
    ymin=id-0.025 , ymax=id+0.025 ), colour="grey50", fill=NA,
    size=0.3)

  }

  # SAVE FINAL OUTPUT 5 - Notes across time and concepts for patients
  pdf(paste0(output_data_folder, "targex_output_patient_time_concept", filename,
   "i_", i, ".pdf"))
  grid.arrange(plot1,plot2, heights=c(2,2), widths=20, ncol=1)
  dev.off()

  print(paste0("Succesfully saved pdf"))
}


################################################################################
###################### (I) OUTPUTS  V  - COLOUR MACRO OUTPUT ###################
################################################################################
### extract notes for the subset of patients for which notes are to be analyised - 
### at most one note per patient
sample <- empi_index_expanded_large

sample_notes <- unique(output[active_note==1], by=c("empi"))[empi %in% 
  sample$empi]

# extract ommited notes
sample_notes <- unique(output[active_note!=1], by=c("empi"))[1:300]

# extract included notes which have never had an encoutner, diagnoses or procedure
sample_notes <- unique(output[active_note==1 & !(empi %in%
unique(c(output_procedures_exp[!is.na(encounter_number) & 
      active_note==1]$empi,output_diagnoses_exp[!is.na(encounter_number)
      & active_note==1]$empi, output_encounters_exp[!is.na(encounter_number) & 
      active_note==1]$empi)) )], by=c("empi"))[1:300]
################################################################################
### extract the relevant full lno notes for the subset of patients
sample_notes_full <- notes.chunks.sent[[1]][0]
for (i in 1:157) {
  temp <- notes.chunks.sent[[i]][noteid %in% sample_notes$note_id]
  sample_notes_full <<- rbind(sample_notes_full,temp)  
}

################################################################################
### merge the targex output with the full text notes
sample_notes_full <- sample_notes_full[order(as.numeric(noteid), 
  as.numeric(sentid))]
sample_notes_full[, sent_count:=.N, by=c("noteid")]
sample_notes_full[, note_text:=paste(sentences,collapse="..."), by=c("noteid")]
sample_notes_full[, note_text:=paste(note_text,"...")]
sample_notes_full <- unique(sample_notes_full, by=c("noteid"))

sample_notes <- sample_notes[order(note_id)]
sample_notes_full <- sample_notes_full[order(noteid)]
sample_notes[, c("note_text","sent_count"):=.(sample_notes_full$note_text,
  sample_notes_full$sent_count)]

sample_notes <- sample_notes[, .SD, .SDcols=c("empi", "lno_date", "note_id",
  "note_text",names(sample_notes)[names(sample_notes) %like% "_sentence"] , 
  "sent_count")]

# SAVE FINAL OUTPUT 6 - CSV file to be used with the colouring macro (formatted 
# according to the colouring macro specifications)
write.csv(sample_notes, file=paste0(output_data_folder, 
  "targex_output_colourcoding", filename, ".csv"), row.names=F)


################################################################################
################################################################################
############################      END      #####################################
################################################################################
################################################################################



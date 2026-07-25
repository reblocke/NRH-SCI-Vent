//NRH SCI Cohort Preprocessing

version 18
clear
set more off
capture log close

args data_dir output_root run_id
if "`data_dir'" == "" local data_dir "Data"
if "`output_root'" == "" local output_root "Results and Figures"

capture confirm file "`data_dir'/Working NRH SCI.csv"
if _rc {
    di as err "Could not find `data_dir'/Working NRH SCI.csv."
    di as err "Place the restricted source CSV in the ignored Data/ directory or pass a data directory as the first argument."
    exit 601
}

if "`run_id'" == "" local results_dir "`output_root'/$S_DATE"
else local results_dir "`output_root'/`run_id'"
local log_dir "`results_dir'/Logs"

capture mkdir "`output_root'"
capture mkdir "`results_dir'" //make new folder for figure output if needed
capture mkdir "`log_dir'" //new folder for stata logs

local nrh_mapping_contract "validation/value_mappings.csv"
local nrh_mapping_diagnostics "`log_dir'/value_mapping_diagnostics.csv"

do "code/lib/data_contracts.do"
do "code/lib/value_mappings.do"
nrh_validate_source_contract using "`data_dir'/Working NRH SCI.csv", ///
    schema("validation/source_schema.csv") ///
    rules("validation/validation_rules.csv") ///
    mode("strict") ///
    log("`log_dir'/source_contract.log")
save "`data_dir'/nrh-sci-raw", replace // 1x command to process the dataset to a stata file

* Data processing
clear
log using "`log_dir'/preprocessing.log", replace text
use "`data_dir'/nrh-sci-raw", clear
capture erase "`nrh_mapping_diagnostics'"

drop pat_mrn_id
drop if missing(age)
assert age >= 0

gen age_decade = age/10
label variable age "Age at inpatient rehabilitation admission"
label variable age_decade "Age (per 10 years)"

rename sex nrh_source_sex
nrh_apply_value_mapping nrh_source_sex, sourcefield("sex") target(sex) ///
    contract("`nrh_mapping_contract'") ///
    diagnostics("`nrh_mapping_diagnostics'")
drop nrh_source_sex
label variable sex "Male sex"

nrh_apply_value_mapping deceased, sourcefield("deceased") target(death) ///
    contract("`nrh_mapping_contract'") ///
    diagnostics("`nrh_mapping_diagnostics'")
drop deceased
label variable death "Died after discharge?"

nrh_apply_value_mapping ribfracturesyorn, ///
    sourcefield("ribfracturesyorn") target(rib_fx) ///
    contract("`nrh_mapping_contract'") ///
    diagnostics("`nrh_mapping_diagnostics'")
drop ribfracturesyorn
label variable rib_fx "Rib Fracture?"

nrh_apply_value_mapping pneumothoraxyorn, ///
    sourcefield("pneumothoraxyorn") target(ptx) ///
    contract("`nrh_mapping_contract'") ///
    diagnostics("`nrh_mapping_diagnostics'")
drop pneumothoraxyorn
label variable ptx "Pneumothorax?"

nrh_apply_value_mapping chesttubeyorn, ///
    sourcefield("chesttubeyorn") target(chest_tube) ///
    contract("`nrh_mapping_contract'") ///
    diagnostics("`nrh_mapping_diagnostics'")
drop chesttubeyorn
label variable chest_tube "Chest Tube Inserted?"

nrh_apply_value_mapping injurylevelreportedpriortorehab, ///
    sourcefield("injurylevelreportedpriortorehab") ///
    target(init_injury_level) contract("`nrh_mapping_contract'") ///
    diagnostics("`nrh_mapping_diagnostics'")
drop injurylevelreportedpriortorehab
label variable init_injury_level "Reported (Prior to Rehab) Injury Level"

// Historical coding note: reason for intubation was not used in the final paper models.
//Cat: GCS low / mental status
//Cat: Respiratory failure (with subclassification for why, if stated)
//Cat: Unknown / on scene

nrh_apply_value_mapping attempttoweanpriortotrachyorn, ///
    sourcefield("attempttoweanpriortotrachyorn") ///
    target(wean_pre_trach) contract("`nrh_mapping_contract'") ///
    diagnostics("`nrh_mapping_diagnostics'")
drop attempttoweanpriortotrachyorn
label variable wean_pre_trach "Attempt to wean before trach?"

nrh_apply_value_mapping attempttoweanoffventpriortotrans, ///
    sourcefield("attempttoweanoffventpriortotrans") ///
    target(wean_pre_trans) contract("`nrh_mapping_contract'") ///
    diagnostics("`nrh_mapping_diagnostics'")
drop attempttoweanoffventpriortotrans
label variable wean_pre_trans "Attempt to wean before transport?"

nrh_apply_value_mapping didpatientdeveloppneumoniapriort, ///
    sourcefield("didpatientdeveloppneumoniapriort") ///
    target(pneumonia_prior) contract("`nrh_mapping_contract'") ///
    diagnostics("`nrh_mapping_diagnostics'")
drop didpatientdeveloppneumoniapriort
label variable pneumonia_prior "Pneumonia prior to t..."

drop column1

nrh_apply_value_mapping asiaclassificationatrehab, ///
    sourcefield("asiaclassificationatrehab") target(level) ///
    contract("`nrh_mapping_contract'") ///
    diagnostics("`nrh_mapping_diagnostics'")
nrh_apply_value_mapping asiaclassificationatrehab, ///
    sourcefield("asiaclassificationatrehab") target(asia_class) ///
    contract("`nrh_mapping_contract'") ///
    diagnostics("`nrh_mapping_diagnostics'")
drop asiaclassificationatrehab
label variable level "Cervical injury level at rehabilitation admission"
label variable asia_class "AIS classification at rehabilitation admission"

generate byte comp_vs_part = .
replace comp_vs_part = 1 if inlist(asia_class, 1, 2)
replace comp_vs_part = 2 if inlist(asia_class, 3, 4)
label define comp_vs_part_label ///
    1 "Complete Motor" 2 "Partial Motor", replace
label values comp_vs_part comp_vs_part_label
label variable comp_vs_part "Complete Motor (AIS A or B) vs Partial Motor (AIS C or AIS D)"

generate byte high_vs_low = .
replace high_vs_low = 1 if inrange(level, 1, 4)
replace high_vs_low = 2 if inrange(level, 5, 8)
label define high_vs_low_label 1 "High" 2 "Low", replace
label values high_vs_low high_vs_low_label
label variable high_vs_low "High (C4 or above) or Low (C5 or below) Cervical Injury?"

generate byte init_high_vs_low = .
replace init_high_vs_low = 1 if inrange(init_injury_level, 1, 4)
replace init_high_vs_low = 2 if inrange(init_injury_level, 5, 8)
label define init_high_vs_low_label 1 "High" 2 "Low", replace
label values init_high_vs_low init_high_vs_low_label
label variable init_high_vs_low "Reported (prior to rehab) as High (C4 or above) or Low (C5 or below) Cervical Injury?"

generate byte reclass_on_arrival = .
replace reclass_on_arrival = 1 if init_injury_level < level & ///
    !missing(init_injury_level, level)
replace reclass_on_arrival = 2 if init_injury_level == level & ///
    !missing(init_injury_level, level)
replace reclass_on_arrival = 3 if init_injury_level > level & ///
    !missing(init_injury_level, level)
label define reclass_on_arrival_label ///
    1 "Down (Lower Level)" 2 "No Change" 3 "Up (Higher Level)", replace
label values reclass_on_arrival reclass_on_arrival_label
label variable reclass_on_arrival "Injury Reclassified at Rehab Admit?"

generate byte level_and_completeness = .
replace level_and_completeness = 0 if high_vs_low == 1 & comp_vs_part == 1
replace level_and_completeness = 1 if high_vs_low == 1 & comp_vs_part == 2
replace level_and_completeness = 2 if high_vs_low == 2 & comp_vs_part == 1
replace level_and_completeness = 3 if high_vs_low == 2 & comp_vs_part == 2
label define level_and_completeness_label ///
    0 "High, Complete Motor" 1 "High, Partial Motor" ///
    2 "Low, Complete Motor" 3 "Low, Partial Motor", replace
label values level_and_completeness level_and_completeness_label
label variable level_and_completeness "Cervical level and motor completeness"

nrh_apply_value_mapping weanoffventduringthedayandcontin, ///
    sourcefield("weanoffventduringthedayandcontin") ///
    target(partial_wean_at_admit) contract("`nrh_mapping_contract'") ///
    diagnostics("`nrh_mapping_diagnostics'")
nrh_apply_value_mapping weanoffventduringthedayandcontin, ///
    sourcefield("weanoffventduringthedayandcontin") ///
    target(wean_during_day) contract("`nrh_mapping_contract'") ///
    diagnostics("`nrh_mapping_diagnostics'")
drop weanoffventduringthedayandcontin
label variable partial_wean_at_admit ///
    "Already weaned during the day before rehab admission?"
label variable wean_during_day "Able to wean during the day?"

foreach nrh_binary in sex death rib_fx ptx chest_tube wean_pre_trach ///
    wean_pre_trans pneumonia_prior partial_wean_at_admit wean_during_day {
    assert missing(`nrh_binary') | inlist(`nrh_binary', 0, 1)
}
assert missing(init_injury_level) | inrange(init_injury_level, 1, 8)
assert missing(level) | inrange(level, 1, 8)
assert missing(asia_class) | inrange(asia_class, 1, 4)
assert missing(comp_vs_part) | inlist(comp_vs_part, 1, 2)
assert missing(high_vs_low) | inlist(high_vs_low, 1, 2)
assert missing(init_high_vs_low) | inlist(init_high_vs_low, 1, 2)
assert missing(reclass_on_arrival) | inrange(reclass_on_arrival, 1, 3)
assert missing(level_and_completeness) | ///
    inrange(level_and_completeness, 0, 3)

replace ifyestopreviousquestionhowmanyda = ///
    ustrupper(ustrtrim(ustrnormalize( ///
        ifyestopreviousquestionhowmanyda, "nfc")))
replace ifyestopreviousquestionhowmanyda = "0" if ///
    ifyestopreviousquestionhowmanyda == "PRIOR TO REHAB"
replace ifyestopreviousquestionhowmanyda = "0" if ///
    ifyestopreviousquestionhowmanyda == "0 (PRIOR TO REHAB)"
replace ifyestopreviousquestionhowmanyda = "0" if ///
    ifyestopreviousquestionhowmanyda == "0- PT WEANED ON ARRIVAL"
destring ifyestopreviousquestionhowmanyda, gen(days_to_daytime_wean)
drop ifyestopreviousquestionhowmanyda
label variable days_to_daytime_wean "Days until weaned during day"

nrh_apply_value_mapping weanoffventforfull24hoursyorn, ///
    sourcefield("weanoffventforfull24hoursyorn") target(wean_24hr) ///
    contract("`nrh_mapping_contract'") ///
    diagnostics("`nrh_mapping_diagnostics'")
drop weanoffventforfull24hoursyorn
label variable wean_24hr "Liberated from IMV? (24hr/d)"

replace v23 = ustrupper(ustrtrim(ustrnormalize(v23, "nfc")))
replace v23 = "1" if v23 == "1 WEANED OFF CPAP AT NIGHT"
destring v23, gen(days_to_24hr_wean) 
drop v23
label variable days_to_24hr_wean "Days until liberated from IMV"

nrh_apply_value_mapping didpatientdecanulateyorn, ///
    sourcefield("didpatientdecanulateyorn") target(decannulate) ///
    contract("`nrh_mapping_contract'") ///
    diagnostics("`nrh_mapping_diagnostics'")
label variable decannulate "Did they decannulate?"
drop didpatientdecanulateyorn

generate byte weaning_outcome = .
replace weaning_outcome = 1 if wean_during_day == 0
replace weaning_outcome = 2 if wean_during_day == 1
replace weaning_outcome = 3 if wean_24hr == 1
replace weaning_outcome = 4 if decannulate == 1
label define weaning_label ///
    1 "Fully Vent Dependent" 2 "Weaned During Day" ///
    3 "Liberated from IMV" 4 "Decannulated", replace
label values weaning_outcome weaning_label
label variable weaning_outcome "Weaning outcome at rehab discharge"

nrh_apply_value_mapping didtheytransfertoicufromrehabyor, ///
    sourcefield("didtheytransfertoicufromrehabyor") ///
    target(rehab_to_icu) contract("`nrh_mapping_contract'") ///
    diagnostics("`nrh_mapping_diagnostics'")
drop didtheytransfertoicufromrehabyor
label variable rehab_to_icu "Did they transfer Rehab->ICU?"

nrh_apply_value_mapping iftheytransferredtoicuwascausere, ///
    sourcefield("iftheytransferredtoicuwascausere") ///
    target(resp_icu_transfer) contract("`nrh_mapping_contract'") ///
    diagnostics("`nrh_mapping_diagnostics'")
drop iftheytransferredtoicuwascausere
label variable resp_icu_transfer "Respiratory cause of ICU transfer?"

nrh_apply_value_mapping didtheydeveloppnumoniaatrehabyor, ///
    sourcefield("didtheydeveloppnumoniaatrehabyor") ///
    target(pna_at_rehab) contract("`nrh_mapping_contract'") ///
    diagnostics("`nrh_mapping_diagnostics'")
drop didtheydeveloppnumoniaatrehabyor
label variable pna_at_rehab "Did they develop pneumonia at rehab?"

tempvar nrh_discharge_normalized
generate strL `nrh_discharge_normalized' = ///
    ustrupper(ustrtrim(ustrnormalize(dischargelocation, "nfc")))
generate byte discharge_to = .
replace discharge_to = 1 if inlist(`nrh_discharge_normalized', ///
    "LTAC", "NEURORESTORITIVE THERAPY")
replace discharge_to = 2 if `nrh_discharge_normalized' == "SNF"
replace discharge_to = 3 if inlist(`nrh_discharge_normalized', ///
    "HH", "HOME HEALTH", "HOME W/ HOME HEALTH")
replace discharge_to = 4 if `nrh_discharge_normalized' == "HOME"
quietly count if !missing(dischargelocation) & missing(discharge_to)
if r(N) {
    di as err "An unlisted discharge value was rejected."
    exit 459
}
label define discharge_label ///
    1 "LTAC" 2 "SNF" 3 "Home w/ Home Health" 4 "Home", replace
label values discharge_to discharge_label
label variable discharge_to "Discharge Location?"
drop dischargelocation

generate byte home_vs_facility = .
replace home_vs_facility = 1 if inlist(discharge_to, 1, 2)
replace home_vs_facility = 2 if inlist(discharge_to, 3, 4)
label variable home_vs_facility "Discharge to Home vs Facility"
label define home_vs_facility_lab ///
    1 "Facility (SNF, LTAC)" 2 "Home (+/- HH)", replace
label values home_vs_facility home_vs_facility_lab

nrh_apply_value_mapping osh, sourcefield("osh") target(outside_hospital) ///
    contract("`nrh_mapping_contract'") ///
    diagnostics("`nrh_mapping_diagnostics'")
label variable outside_hospital "Admit from Outside Hospital?"
drop osh

foreach nrh_binary in wean_24hr decannulate rehab_to_icu ///
    resp_icu_transfer pna_at_rehab outside_hospital {
    assert missing(`nrh_binary') | inlist(`nrh_binary', 0, 1)
}
assert missing(weaning_outcome) | inrange(weaning_outcome, 1, 4)
assert missing(discharge_to) | inrange(discharge_to, 1, 4)
assert missing(home_vs_facility) | inlist(home_vs_facility, 1, 2)

gen injury_date = date(strtrim(dateofinjury), "MD20Y")
format injury_date %td
drop dateofinjury
label variable injury_date "Date of Injury"

//cause of death
//daysfromdischargetodeath

gen date_death = date(strtrim(tod), "MD20Y")
label variable date_death "Date of Death"
drop tod
format date_death %td

gen date_discharge = injury_date + daysfrominjurytoadmissiontorehab + daystodischargefromrehab
format date_discharge %td
label variable date_discharge "Date of Discharge"

gen date_admit = injury_date + daysfrominjurytoadmissiontorehab
format date_admit %td
label variable date_admit "Date of Admit"

local nrh_admin_censor_date_iso "2023-09-30"
local nrh_administrative_censor_date = ///
    daily("`nrh_admin_censor_date_iso'", "YMD")
if missing(`nrh_administrative_censor_date') {
    di as err "The version-controlled administrative censor date is invalid."
    exit 459
}

gen time_to_censor_death = date_death - date_admit if death == 1 & ///
    !missing(date_death, date_admit)
replace time_to_censor_death = ///
    `nrh_administrative_censor_date' - date_admit if death == 0 & ///
    !missing(date_admit)
label variable time_to_censor_death "Days of follow-up (or until death)"

gen time_to_censor_death_dc = date_death - date_discharge if death == 1 & ///
    !missing(date_death, date_discharge)
replace time_to_censor_death_dc = ///
    `nrh_administrative_censor_date' - date_discharge if death == 0 & ///
    !missing(date_discharge)
assert missing(time_to_censor_death) | time_to_censor_death >= 0
assert missing(time_to_censor_death_dc) | time_to_censor_death_dc >= 0

save "`data_dir'/nrh-sci-cleaned", replace
log close

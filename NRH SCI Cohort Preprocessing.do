//NRH SCI Cohort Preprocessing

clear
cd "/Users/blocke/Box Sync/Residency Personal Files/Scholarly Work/Locke Research Projects/NRH SCI Data/" 

program define datetime 
end

capture mkdir "Results and Figures"
capture mkdir "Results and Figures/$S_DATE/" //make new folder for figure output if needed
capture mkdir "Results and Figures/$S_DATE/Logs/" //new folder for stata logs


import delimited using "Data/Working NRH SCI.csv", colrange(1:36) //, sheet("All") firstrow //case(lower) 
save nrh-sci-raw, replace // 1x command to process the dataset to a stata file

* Data processing
clear
capture log close
log using "Results and Figures/$S_DATE/Logs/temp.log", append
use nrh-sci-raw

drop pat_mrn_id
drop if missing(age)
label define binary_label 0 "N" 1 "Y"

gen age_decade = age/10
label variable age_decade "Age (per 10 years)"

replace sex = "M" if sex == "male"
encode sex, gen(male) label(binary_label)
drop sex
rename male sex
label variable sex "Sex"

encode deceased, gen(death) label(binary_label)
drop deceased
label variable death "Died after discharge?"

encode ribfracturesyorn, gen(rib_fx) label(binary_label)
drop ribfracturesyorn
label variable rib_fx "Rib Fracture?"

encode pneumothoraxyorn, gen(ptx) label(binary_label)
drop pneumothoraxyorn
label variable ptx "Pneumothorax?"

replace chesttubeyorn = "Y" if chesttubeyorn != "N"
encode chesttubeyorn, gen(chest_tube) label(binary_label)
drop chesttubeyorn
label variable chest_tube "Chest Tube Inserted?"

label define level_cat 1 "C1" 2 "C2" 3 "C3" 4 "C4" 5 "C5" 6 "C6" 7 "C7"

replace injurylevelreportedpriortorehab = strtrim(strupper(injurylevelreportedpriortorehab))
encode injurylevelreportedpriortorehab, gen(init_injury_level) label(level_cat)
drop injurylevelreportedpriortorehab
label variable init_injury_level "Reported (Prior to Rehab) Injury Level"

//TODO: coding of reason for intubation
//Cat: GCS low / mental status
//Cat: Respiratory failure (with subclassification for why, if stated)
//Cat: Unknown / on scene

replace attempttoweanpriortotrachyorn = strtrim(strupper(attempttoweanpriortotrachyorn))
replace attempttoweanpriortotrachyorn = "N" if attempttoweanpriortotrachyorn == "NO"
replace attempttoweanpriortotrachyorn = "Y" if attempttoweanpriortotrachyorn == "Y - EXTUBATED ON 8/20  BUT PT STRIDORED AND C.A."
encode attempttoweanpriortotrachyorn, label(binary_label) gen(wean_pre_trach)
replace wean_pre_trach = . if !inlist(wean_pre_trach, 0,1)
drop attempttoweanpriortotrachyorn
label variable wean_pre_trach "Attempt to wean before trach?"

replace attempttoweanoffventpriortotrans = strtrim(strupper(attempttoweanoffventpriortotrans))
replace attempttoweanoffventpriortotrans = "Y" if attempttoweanoffventpriortotrans == "Y (ABLE TO WEANE TO VENT AT NIGHT)"
encode attempttoweanoffventpriortotrans , label(binary_label) gen(wean_pre_trans)
replace wean_pre_trans = . if !inlist(wean_pre_trans, 0,1)
drop attempttoweanoffventpriortotrans 
label variable wean_pre_trans "Attempt to wean before transport?"

replace didpatientdeveloppneumoniapriort = strtrim(strupper(didpatientdeveloppneumoniapriort))
replace didpatientdeveloppneumoniapriort = "Y" if didpatientdeveloppneumoniapriort == "YES"
encode didpatientdeveloppneumoniapriort, label(binary_label) gen(pneumonia_prior)
drop didpatientdeveloppneumoniapriort
label variable pneumonia_prior "Pneumonia prior to t..."

drop column1

//split apart variable
replace asiaclassificationatrehab = strtrim(strupper(asiaclassificationatrehab))
split asiaclassificationatrehab, gen(level_class)
replace level_class2 = "AIS" if level_class2 == "AIA"
encode level_class1, generate(level) label(level_cat)
label variable level "Level of Injury"
drop level_class1
generate asia_class_temp = level_class2 + " " + level_class3
replace asia_class_temp = "AIS D" if asia_class_temp == "CENTRAL CORD"
encode asia_class_temp, gen(asia_class)
//list level asia_class level_class2 level_class3 asia_class_temp asiaclassificationatrehab
drop level_class2 level_class3 asia_class_temp asiaclassificationatrehab //keep level and asia_class
label variable asia_class "ASIA classification"

generate comp_v_part = ""
replace comp_v_part = "Complete Motor" if asia_class == 1 //AIS A
replace comp_v_part = "Complete Motor" if asia_class == 2 //AIS B
replace comp_v_part = "Partial Motor" if asia_class == 3 //AIS C
replace comp_v_part = "Partial Motor" if asia_class == 4 //AIS 4  aka Central Cord
encode comp_v_part, gen(comp_vs_part)  //1 = comp, 2 = partial
drop comp_v_part
label variable comp_vs_part "Complete Motor (AIS A or B) vs Partial Motor (AIS C or AIS D)"

generate high_v_low = ""
replace high_v_low = "High" if level == 1
replace high_v_low = "High" if level == 2
replace high_v_low = "High" if level == 3
replace high_v_low = "High" if level == 4
replace high_v_low = "Low" if level == 5
replace high_v_low = "Low" if level == 6
replace high_v_low = "Low" if level == 7
encode high_v_low, gen(high_vs_low)  //High = 1, Low = 2
drop high_v_low
label variable high_vs_low "High (C4 or above) or Low (C5 or below) Cervical Injury?"

generate init_high_v_low = ""
replace init_high_v_low = "High" if init_injury_level == 1
replace init_high_v_low = "High" if init_injury_level == 2
replace init_high_v_low = "High" if init_injury_level == 3
replace init_high_v_low = "High" if init_injury_level == 4
replace init_high_v_low = "Low" if init_injury_level == 5
replace init_high_v_low = "Low" if init_injury_level == 6
replace init_high_v_low = "Low" if init_injury_level == 7
encode init_high_v_low, gen(init_high_vs_low)
drop init_high_v_low
label variable init_high_vs_low "Reported (prior to rehab) as High (C4 or above) or Low (C5 or below) Cervical Injury?"

generate reclass_arrival = ""
replace reclass_arrival = "Up (Higher Level)" if init_injury_level > (level) 
replace reclass_arrival = "No Change" if init_injury_level == (level)
replace reclass_arrival = "Down (Lower Level)" if init_injury_level < (level)
encode reclass_arrival, gen(reclass_on_arrival)
label variable reclass_on_arrival "Injury Reclassified at Rehab Admit?"
drop reclass_arrival

gen level_completeness = ""
label define level_comp_label 0 "High, Complete Injury" 1 "High, Incomplete Injury" 2 "Low, Complete Injury" 3 "Low, Incomplete Injury"
replace level_completeness = "High, Complete Injury" if (high_vs_low == 1) & (comp_vs_part == 1)
replace level_completeness = "High, Incomplete Injury" if (high_vs_low == 1) & (comp_vs_part == 2)
replace level_completeness = "Low, Complete Injury" if (high_vs_low == 2) & (comp_vs_part == 1)
replace level_completeness = "Low, Incomplete Injury" if (high_vs_low == 2) & (comp_vs_part == 2)
encode level_completeness, gen(level_and_completeness) label(level_comp_label)
drop level_completeness

gen partial_wean_admit = "N"
replace weanoffventduringthedayandcontin = strtrim(strupper(weanoffventduringthedayandcontin))
replace partial_wean_admit = "Y" if weanoffventduringthedayandcontin == "Y (PRIOR TO REHAB)"
replace weanoffventduringthedayandcontin = "Y" if weanoffventduringthedayandcontin == "Y (PRIOR TO REHAB)"
replace partial_wean_admit = "Y" if weanoffventduringthedayandcontin == "Y PRIOR TO ADMISSION"
replace weanoffventduringthedayandcontin = "Y" if weanoffventduringthedayandcontin == "Y PRIOR TO ADMISSION"
encode partial_wean_admit, generate(partial_wean_at_admit) label(binary_label)
label variable partial_wean_at_admit "Already Weaned During the Day at Admit?"
drop partial_wean_admit
replace weanoffventduringthedayandcontin = "Y" if weanoffventduringthedayandcontin == "Y."
encode weanoffventduringthedayandcontin, gen(wean_during_day) label(binary_label)
drop weanoffventduringthedayandcontin
label variable wean_during_day "Able to wean during the day?"

replace ifyestopreviousquestionhowmanyda = strtrim(ifyestopreviousquestionhowmanyda)
replace ifyestopreviousquestionhowmanyda = "0" if ifyestopreviousquestionhowmanyda == "Prior to rehab"
replace ifyestopreviousquestionhowmanyda = "0" if ifyestopreviousquestionhowmanyda == "0 (Prior to rehab)"
replace ifyestopreviousquestionhowmanyda = "0" if ifyestopreviousquestionhowmanyda == "0- pt weaned on arrival"
destring ifyestopreviousquestionhowmanyda, gen(days_to_daytime_wean)
drop ifyestopreviousquestionhowmanyda
label variable days_to_daytime_wean "Days until weaned during day"

replace weanoffventforfull24hoursyorn = strtrim(strupper(weanoffventforfull24hoursyorn))
replace weanoffventforfull24hoursyorn = "N" if weanoffventforfull24hoursyorn == "N - BIPAP AT NIGHT"
encode weanoffventforfull24hoursyorn, gen(wean_24hr) label(binary_label)
drop weanoffventforfull24hoursyorn
label variable wean_24hr "Liberated from IMV? (24hr/d)"

replace v23 = strtrim(v23)
replace v23 = "1" if v23 == "1 weaned off CPAP at night"
destring v23, gen(days_to_24hr_wean) 
drop v23
label variable days_to_24hr_wean "Days until liberated from IMV"

replace didpatientdecanulateyorn = strtrim(strupper(didpatientdecanulateyorn))
replace didpatientdecanulateyorn = "N" if didpatientdecanulateyorn == "N."
encode didpatientdecanulateyorn, gen(decannulate) label(binary_label)
label variable decannulate "Did they decannulate?"
drop didpatientdecanulateyorn

gen weaning_outcome_str = ""
replace weaning_outcome_str = "Fully Vent Dependent" if wean_during_day == 0 
replace weaning_outcome_str = "Weaned During Day" if wean_during_day == 1
replace weaning_outcome_str = "Liberated from IMV" if wean_24hr == 1
replace weaning_outcome_str = "Decannulated" if decannulate == 1
label define weaning_label 1 "Fully Vent Dependent" 2 "Weaned During Day" 3 "Liberated from IMV" 4 "Decannulated"
encode weaning_outcome_str, label(weaning_label) gen(weaning_outcome)
drop weaning_outcome_str
label variable weaning_outcome "Weaning outcome at rehab discharge"

replace didtheytransfertoicufromrehabyor = strtrim(strupper(didtheytransfertoicufromrehabyor))
replace didtheytransfertoicufromrehabyor = "Y" if didtheytransfertoicufromrehabyor == "Y  - SICU"
encode didtheytransfertoicufromrehabyor, gen(rehab_to_icu) label(binary_label)
drop didtheytransfertoicufromrehabyor
label variable rehab_to_icu "Did they transfer Rehab->ICU?"

replace iftheytransferredtoicuwascausere = strtrim(strupper(iftheytransferredtoicuwascausere))
replace iftheytransferredtoicuwascausere = "N" if iftheytransferredtoicuwascausere == "N (DIAPHRAM PACER PLACEMENT)"
encode iftheytransferredtoicuwascausere, gen(resp_icu_transfer) label(binary_label)
drop iftheytransferredtoicuwascausere
label variable resp_icu_transfer "Respiratory cause of ICU transfer?"

replace didtheydeveloppnumoniaatrehabyor = strtrim(strupper(didtheydeveloppnumoniaatrehabyor))
encode didtheydeveloppnumoniaatrehabyor, gen(pna_at_rehab) label(binary_label)
drop didtheydeveloppnumoniaatrehabyor
label variable pna_at_rehab "Did they develop pneumonia at rehab?"

label define discharge_label 4 "Home" 3 "Home w/ Home Health" 2 "SNF" 1 "LTAC"
replace dischargelocation = strtrim(strupper(dischargelocation))
replace dischargelocation = "HOME W/ HOME HEALTH" if dischargelocation == "HOME HEALTH"
replace dischargelocation = "HOME W/ HOME HEALTH" if dischargelocation == "HH"
replace dischargelocation = "LTAC" if dischargelocation == "NEURORESTORITIVE THERAPY"
replace dischargelocation = "Home" if dischargelocation == "HOME"
replace dischargelocation = "Home w/ Home Health" if dischargelocation == "HOME W/ HOME HEALTH"

encode dischargelocation, label(discharge_label) gen(discharge_to)
label variable discharge_to "Discharge Location?"
drop dischargelocation

recode discharge_to (1=1) (2=1) (3=2) (4=2), generate(home_vs_facility)
label variable home_vs_facility "Discharge to Home vs Facility"
label define home_vs_facility_lab 1 "Facility (SNF, LTAC)" 2 "Home (+/- HH)"
label values home_vs_facility home_vs_facility_lab

replace osh  = strtrim(strupper(osh))
encode osh, gen(outside_hospital) label(binary_label)
label variable outside_hospital "Admit from Outside Hospital?"
drop osh

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

/* Old? */ 

gen time_to_censor_death = date_death - date_admit if death == 1
replace time_to_censor_death = date("09/30/2023", "MDY") - date_admit if death == 0
label variable time_to_censor_death "Days of follow-up (or until death)"

gen time_to_censor_death_dc = date_death - date_discharge if death == 1
replace time_to_censor_death_dc = date("09/30/2023", "MDY") - date_discharge if death == 0

save nrh-sci-cleaned, replace

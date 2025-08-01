* NRH SCI Cohort Analysis

clear
cd "/Users/blocke/Box Sync/Residency Personal Files/Scholarly Work/Locke Research Projects/NRH SCI Data" //Mac

program define datetime 
end

capture mkdir "Results and Figures"
capture mkdir "Results and Figures/$S_DATE/" //make new folder for figure output if needed
capture mkdir "Results and Figures/$S_DATE/Logs/" //new folder for stata logs

* Data Analysis

clear
capture log close
log using "Results and Figures/$S_DATE/Logs/temp.log", append

clear
use nrh-sci-cleaned

/* Analysis: 
Exposure: Age, Injury (high vs low), Completeness (high vs low) 

Outcome 1: Weaning Outcome Category [raw; then ordinal regression; then Fine-Gray]
Outcome 2: Discharge Location [raw; then ordinal regression]
Outcome 3: Death (separate info) [raw; then K-M]
*/

/* Exclude if partial wean at admit, n=3 */ 
drop if partial_wean_at_admit == 1

/* Figure 1 manually made by CF */ 
/* Table 1  - Demographics & Pre-Rehab Course */ 

table1_mc, ///
vars( ///
sex cat %4.0f \ ///
age conts %4.0f \ ///
level cat %4.0f \ ///
asia_class cat %4.0f \ ///
reclass_on_arrival cat %4.0f \ ///
outside_hospital bin %4.0f \ ///
daysfrominjurytointubation conts %4.0f \ ///
daysfromintubationtotrach conts %4.0f \ ///
daysfrominjurytoadmissiontorehab conts %4.0f \ ///
) ///
percent_n percsign("%") iqrmiddle(",") sdleft(" (±") sdright(")") onecol total(before) ///
saving("Results and Figures/$S_DATE/Table 1 - PreRehab.xlsx", replace)

//Overall weaning course
table1_mc, ///
vars( ///
daystodischargefromrehab conts %4.0f \ ///
wean_during_day bin %4.0f \ ///
days_to_daytime_wean conts %4.0f \ ///
wean_24hr bin %4.0f \ ///
days_to_24hr_wean conts %4.0f \ ///
decannulate bin %4.0f \ ///
daysfromadmissiontorehabtodecanu conts %4.0f \ ///
discharge_to cat %4.0f \ ///
) ///
percent_n percsign("%") iqrmiddle(",") sdleft(" (±") sdright(")") onecol total(before) ///
saving("Results and Figures/$S_DATE/Table aux - Rehab milestones overall.xlsx", replace)

/* 
Table: Discharge Location, Weaning Status by Injury 
*/

table1_mc, by(high_vs_low) ///
vars( ///
weaning_outcome cat %4.0f \ ///
daystodischargefromrehab conts %4.0f \ ///
discharge_to cat %4.0f \ ///
death bin %4.0f \ ///
) ///
percent_n percsign("%") iqrmiddle(",") sdleft(" (±") sdright(")") onecol ///
saving("Results and Figures/$S_DATE/Table 2a - Outcomes by Level.xlsx", replace)

/* 
Table: Discharge Locations by weaning cat
*/
table1_mc, by(weaning_outcome) ///
vars( ///
daystodischargefromrehab conts %4.0f \ ///
discharge_to cat %4.0f \ ///
death bin %4.0f \ ///
) ///
total(before) percent_n percsign("%") iqrmiddle(",") sdleft(" (±") sdright(")") onecol ///
saving("Results and Figures/$S_DATE/Table 2b - Outcomes by Weaning Cat.xlsx", replace)


/*****************************************************************
Figure 2- PATIENT-DAY PANEL 
*****************************************************************/
preserve
gen id   = _n
local nd = 130
expand `nd'
bys id : gen day = _n - 1

/* 1) set all in-hospital states first  */
/* 1) in-hospital states (0–3) */
gen byte state = 7  
replace state = 6 if day>=daystodischargefromrehab & weaning_outcome==1  // vent-dep

replace state = 5 if day>=days_to_daytime_wean & day<days_to_24hr_wean
replace state = 4 if day>=daystodischargefromrehab & weaning_outcome==2  // daytime

replace state = 3 if day>=days_to_24hr_wean  & day<daysfromadmissiontorehabtodecanu
replace state = 2 if day>=daystodischargefromrehab & weaning_outcome==3  // 24 h off

replace state = 1 if day>=daysfromadmissiontorehabtodecanu & day<daystodischargefromrehab
replace state = 0 if day>=daystodischargefromrehab & weaning_outcome==4  // decannulated

/* 3) relabel in the correct numeric order */
label define statelab ///
  7 "Fully Ventilator Dependent"   ///
  6 "Discharged Fully Vent. Dependent" ///
  5 "Daytime Wean" ///
  4 "Discharged Daytime Wean"  ///  
  3 "Complete Wean"    ///
  2 "Discharged Complete Wean"    ///  
  1 "Decannulated" ///
  0 "Discharged Decannulated", replace
label values state statelab

/* 4) now draw your stackedcount */
stackedcount state day, ///
    xlabel(0(10)130)                   ///
    ytitle("Weaning & Discharge Status")        ///
    xtitle("Day from Rehab. Admission")        ///
	legend(rows(8) position(9) size(small)) ///
    scheme(white_w3d)

graph export "Results and Figures/$S_DATE/Fig 2 - stacked_states.png", as(png) replace ///
    width(2700) height(1500)
restore



/* Supplemental Figure - Time to weaning outcomes */ 
/*
Subdistribution Hazards (Fine-Gray) 

Figure 2 a-c: 
Survival analysis / K-M stratified by low vs high, complete vs incomplete 

A: survival rehab-day wean
B: survival rehab-complete wean
C: survival rehab to decannulation

Note: discharge is a competing risk, thus it might be appropriate to investigate using a fine-gray model . 

There were no inpatient-rehab deaths to be counted as competing events. 
*/ 

/* "Survival" Until Daytime Wean */

preserve 

// 1 = weaned during day
// 2 = discharged before weaning (competing event)
// All other cases (e.g., still in rehab and not weaned): remain missing (censored)
gen event_type = .
replace event_type = 1 if wean_during_day == 1
replace event_type = 2 if daystodischargefromrehab < days_to_daytime_wean & missing(event_type)

gen time_to_event = min(days_to_daytime_wean, daystodischargefromrehab)
stset time_to_event, failure(event_type == 1)
stcrreg high_vs_low c.age_decade i.comp_vs_part, compete(event_type == 2)

stcurve, cif at1(high_vs_low=2) at2(high_vs_low=1) ///
lwidth(thick thick)                                             ///
	legend(order(1 "Low Cervical (C5 & Below)" 2 "High Cervical (C4 & Above)") position(4) ring(0) rows(2) size(medsmall)) ///
	title("Daytime Weaning", size(large)) ///
	xtitle("Day of Rehab Stay", size(medlarge)) ///
	ytitle("Cumulative Incidence of Daytime Wean", size(medlarge)) ///
	xlabel(0(10)80, labsize(medlarge)) ///
	ylabel(0(.1)1, labsize(medlarge)) ///
	range(0 80) 
graph export "Results and Figures/$S_DATE/CIF_DayWean_Level.png", as(png) replace
graph save "CIF_day_wean.gph", replace
restore

/* "Survival" Until 24H Wean */
preserve

// 1 = weaned within 24 h
// 2 = discharged before 24 h wean (competing event)
// All other cases (e.g., still in rehab and not weaned): remain missing (censored)
gen event_type = .
replace event_type = 1 if wean_24hr == 1
replace event_type = 2 if daystodischargefromrehab < days_to_24hr_wean & missing(event_type)

// analysis-time variable = first of the two possible times
gen time_to_event = min(days_to_24hr_wean, daystodischargefromrehab)

// declare survival data and fit competing-risk model
stset time_to_event, failure(event_type == 1)
stcrreg high_vs_low c.age_decade i.comp_vs_part, compete(event_type == 2)

stcurve, cif at1(high_vs_low=2) at2(high_vs_low=1)                  ///
        lwidth(thick thick)                                         ///
        legend(order(1 "Low Cervical (C5 & Below)"                  ///
                     2 "High Cervical (C4 & Above)") position(10)    ///
               ring(0) rows(2) size(medsmall))                      ///
        title("24-Hour Weaning", size(large)) ///
        xtitle("Day of Rehab Stay",   size(medlarge))               ///
        ytitle("Cumulative Incidence of 24-h Wean", size(medlarge)) ///
        xlabel(0(10)80, labsize(medlarge))                          ///
        ylabel(0(.1)1, labsize(medlarge))                           ///
        range(0 80)    
graph export "Results and Figures/${S_DATE}/CIF_24hrWean_Level.png", ///
            as(png) replace
graph save "CIF_24_wean.gph", replace
restore


/* "Survival" Until Decannulation */ 
preserve

// 1 = decannulated
// 2 = discharged before decannulation (competing event)
// All other cases remain missing (censored)
gen event_type = .
replace event_type = 1 if decannulate == 1
replace event_type = 2 if daystodischargefromrehab < daysfromadmissiontorehabtodecanu & missing(event_type)

// analysis-time variable = sooner of decannulation or discharge
gen time_to_event = min(daysfromadmissiontorehabtodecanu, daystodischargefromrehab)

// declare survival data and fit competing-risk model
stset  time_to_event, failure(event_type == 1)
stcrreg high_vs_low c.age_decade i.comp_vs_part, compete(event_type == 2)

// cumulative-incidence plot, same styling
stcurve, cif at1(high_vs_low=2) at2(high_vs_low=1)                  ///
        lwidth(thick thick)                                         ///
        legend(order(1 "Low Cervical (C5 & Below)"                  ///
                     2 "High Cervical (C4 & Above)") position(10)    ///
               ring(0) rows(2) size(medsmall))                      ///
        title("Tracheostomy Decannulation", size(large)) ///
        xtitle("Day of Rehab Stay", size(medlarge))                 ///
        ytitle("Cumulative Incidence of Decannulation", size(medlarge)) ///
        xlabel(0(10)80, labsize(medlarge))                          ///
        ylabel(0(.1)1, labsize(medlarge))                           ///
		text(0.7 0 "*One patient with a high CSCI decannulated at day 179", placement(e) size(small)) ///
        range(0 80)

// export figure
graph export "Results and Figures/${S_DATE}/CIF_Decannulation_Level.png", ///
            as(png) replace
graph save "CIF_decannulation.gph", replace
restore

graph combine CIF_day_wean.gph CIF_24_wean.gph CIF_decannulation.gph, ///
	cols(1) /// 
	xsize(5) ysize(10)
graph export "Results and Figures/$S_DATE/Supp Figure - CIFs for milestones.png", name("Graph") replace



/*
Figure 3: 
Coefplot of regression model: ordinal logit [level of weaning] with predictors: high_v_low, comp_v_inc, c.age
*/ 
ologit weaning_outcome c.age_decade i.comp_vs_part i.high_vs_low, or
estimates store ord_weaning_outcome_reg

ologit discharge_to c.age_decade i.comp_vs_part i.high_vs_low, or
estimates store ord_discharge_to_reg

coefplot ord_weaning_outcome_reg, bylabel("Weaning Outcome") || ord_discharge_to_reg, bylabel("Discharge Location") ||, eform xscale(log) xline(1) xlabel(0.25 0.5 1 2 4 8 16) xscale(extend) xtitle("Odds Ratio of a better weaning or discharge category" , size(small)) yscale(extend) ciopts(recast(rcap) lwidth(thick)) mlabel(string(@b,"%9.2f") + " [ " + string(@ll,"%9.2f") + " - " + string(@ul,"%9.2f") + " ] " + cond(@pval<.001, "***", cond(@pval<.01, "**", cond(@pval<.05, "*", "")))) mlabsize(medsmall) mlabposition(12) mlabgap(*1) headings(age_decade = "{bf:Age}" 2.comp_vs_part = "{bf:Injury} (vs complete)" 2.high_vs_low = "{bf:Level} (vs High)") scheme(white_tableau) text(3 0.5 "Favors Worse" "Category" 3 5.0 "Favors Better" "Category", size(small) color(gs9))
graph export "Results and Figures/$S_DATE/Fig 3 - Ordinal Regressions.png", as(png) name("Graph") replace


/* Supplemental Figure - does weaning status influence discharge location? */ 

ologit discharge_to ib1.weaning_outcome c.age, or
estimates store age_adj_wean_to_discharge
coefplot age_adj_wean_to_discharge, eform xscale(log) xline(1) xlabel(0.061 0.125 0.25 0.5 1 2 4 8 16 32 64 128) xscale(extend) xtitle("Odds Ratio of a lower level of care (LOC) discharge location" , size(small)) yscale(extend) ciopts(recast(rcap) lwidth(thick)) mlabel(string(@b,"%9.2f") + " [ " + string(@ll,"%9.2f") + " - " + string(@ul,"%9.2f") + " ] " + cond(@pval<.001, "***", cond(@pval<.01, "**", cond(@pval<.05, "*", "")))) mlabsize(medsmall) mlabposition(12) mlabgap(*1) scheme(white_tableau) text(5 0.11 "Favors Higher" "LOC Discharge" 5 9.9 "Favors Lower" "LOC Discharge", size(small) color(gs9)) headings(age = "{bf:Age} (per add'n year)" 2.weaning_outcome = "{bf:Weaning Outcome?} (vs 24h Vent)")
graph export "Results and Figures/$S_DATE/Supplement - Dispo by Weaning Status  Regression.png", as(png) name("Graph") replace



/* 
Supplemental Table : age, injury, and weaning outcome by death status
*/ 
table1_mc, by(death) ///
vars( ///
weaning_outcome cat %4.0f \ ///
discharge_to cat %4.0f \ ///
time_to_censor_death conts %4.0f \ ///
time_to_censor_death_dc conts %4.0f \ ///
) ///
percent_n percsign("%") iqrmiddle(",") sdleft(" (±") sdright(")") onecol ///
saving("Results and Figures/$S_DATE/Supplement - Outcome and discharge by Death Status.xlsx", replace)


//Figure 4 Mortality outcomes

/* Survival by Discharge Location */ 
stset time_to_censor_death, failure(death==1)
sts test discharge_to, logrank 

sts graph, by(discharge_to) tmax(2190) ///
 plotopts(lwidth(thick)) ///
 risktable(0(365)2190, order(1 "LTAC" 2 "SNF" 3 "Home with HH" 4 "Home") size(medlarge) title("Number Eligible", size(medlarge))) ///
 text(0.2 30 "Log-Rank {it:p} = 0.04", placement(e) size(medlarge)) ///
 xlabel(0(365)2190, labsize(medlarge)) ///
 ylabel(, labsize(medlarge)) ///
 xtitle("Day From Rehab Admission", size(medlarge)) ///
 ytitle("Survival", size(medlarge)) ///
 legend(order(1 "LTAC" 2 "SNF" 3 "Home w HH" 4 "Home") position(4) ring(0) rows(2) size(medsmall)) ///
 title("Mortality by Discharge Location", size(large))
graph export "Results and Figures/$S_DATE/Supp - KM Death by Dispo.png", as(png) name("Graph") replace
graph save "KM_death_by_dispo.gph", replace
stcox ib1.discharge_to 


stset time_to_censor_death, failure(death==1)
sts test weaning_outcome, logrank 

sts graph, by(weaning_outcome) tmax(2190) ///
 plotopts(lwidth(thick)) ///
 risktable(0(365)2190, order(1 "Full Vent Support" 2 "Noct. Vent Support" 3 "No Vent Support" 4 "Decannulated") size(medlarge) title("Number Eligible", size(medlarge))) ///
 text(0.2 30 "Log-Rank {it:p} = 0.37", placement(e) size(medlarge)) ///
 xlabel(0(365)2190, labsize(medlarge)) ///
 ylabel(, labsize(medlarge)) ///
 xtitle("Day From Rehab Admission", size(medlarge)) ///
 ytitle("Survival", size(medlarge)) ///
 legend(order(1 "Full Vent" 2 "Noct Vent" 3 "No Vent" 4 "Decannulated") position(4) ring(0) rows(2) size(medsmall)) ///
 title("Mortality by Discharge Ventilator Support", size(large))
graph export "Results and Figures/$S_DATE/Supp - KM Death by Wean.png", as(png) name("Graph") replace
graph save "KM_death_by_wean.gph", replace
stcox ib1.weaning_outcome

graph combine KM_death_by_dispo.gph KM_death_by_wean.gph, ///
	cols(2) /// 
	xsize(10) ysize(5)
graph export "Results and Figures/$S_DATE/Figure 4 - KMs for death.png", name("Graph") replace

log close


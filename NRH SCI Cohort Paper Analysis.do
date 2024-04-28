* NRH SCI Cohort Analysis
* 2023-3-30; Update 9/21/23; Update 4/27/24 Brian Locke

clear
cd "/Users/blocke/Box Sync/Residency Personal Files/Scholarly Work/Locke Research Projects/NRH SCI Data" //Mac
//cd "C:\Users\reblo\Box\Residency Personal Files\Scholarly Work\Locke Research Projects\NRH SCI Data" //PC

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

Outcome 1: Weaning Outcome Category [raw; then ordinal regression; then K-M]
Outcome 2: Discharge Location [raw; then ordinal regression]
Outcome 3: Death (separate info) [raw; then K-M]
*/


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
partial_wean_at_admit cat %4.0f \ ///
) ///
percent_n percsign("%") iqrmiddle(",") sdleft(" (±") sdright(")") onecol total(before) ///
saving("Results and Figures/$S_DATE/Table 1 - PreRehab.xlsx", replace)


/* 
Table 2a: Discharge Location by Injury 
*/

table1_mc, by(level_and_completeness) ///
vars( ///
daystodischargefromrehab conts %4.0f \ ///
discharge_to cat %4.0f \ ///
death bin %4.0f \ ///
) ///
percent_n percsign("%") iqrmiddle(",") sdleft(" (±") sdright(")") onecol ///
saving("Results and Figures/$S_DATE/Table 2a - Outcomes by Weaning Cat.xlsx", replace)

/* 
Table 2b: Discharge Locations by weaning cat
*/
table1_mc, by(weaning_outcome) ///
vars( ///
daystodischargefromrehab conts %4.0f \ ///
discharge_to cat %4.0f \ ///
death bin %4.0f \ ///
) ///
percent_n percsign("%") iqrmiddle(",") sdleft(" (±") sdright(")") onecol ///
saving("Results and Figures/$S_DATE/Table 2b - Outcomes by Weaning Cat.xlsx", replace)


/* Figure 2 - Time to weaning outcomes */ 
/*
Figure 2 a-c: 
Survival analysis / K-M stratified by low vs high, complete vs incomplete 

A: survival rehab-day wean
B: survival rehab-complete wean
C: survival rehab to decannulation


Note: discharge is a competing risk, thus it might be appropriate to investigate using a fine-gray model . 

*/ 

/* "Survival" Until Daytime Wean */
preserve
drop if partial_wean_at_admit == 2 //exclude those already partially weaned.
tab high_vs_low wean_during_day
hist days_to_daytime_wean, by(high_vs_low)
gen censored_days_to_daytime_wean = days_to_daytime_wean

/* change comment for sensitivity analysis with noninformative censoring (discharge): */ 
replace censored_days_to_daytime_wean = 180 if missing(days_to_daytime_wean)
//replace sens_cens_days_to_daytime_wean = daystodischargefromrehab  if missing(days_to_daytime_wean)
//note("Sensitivity analysis: Non-informative Censoring")

stset censored_days_to_daytime_wean, failure(wean_during_day==1)
sts test high_vs_low, logrank 
sts graph, by(high_vs_low) failure tmax(90) ///
 plotopts(lwidth(thick)) ///
 risktable(0(10)90, order(2 "Low Cervical" 1 "High Cervical") size(medlarge) title("Number Eligible", size(medlarge))) ///
 text(0.9 5 "Log-Rank {it:p} = 0.005", placement(e) size(medlarge)) ///
 xlabel(0(10)90, labsize(medlarge)) ///
 ylabel(, labsize(medlarge)) ///
 xtitle("Day of Rehab Stay", size(medlarge)) ///
 ytitle("Portion Weaned During Day", size(medlarge)) ///
 legend(order(1 "High Cervical (C4 & Above)" 2 "Low Cervical (C5 & Below)") position(4) ring(0) rows(2) size(medsmall)) ///
 title("When Daytime Wean Occured", size(large))
graph export "Results and Figures/$S_DATE/Unadjusted KM Day-wean Level.png", as(png) name("Graph") replace
graph save "KM_day_wean.gph", replace
stcox high_vs_low
restore

/* "Survival" Until 24H Wean */
preserve
drop if partial_wean_at_admit == 2 //exclude those already partially weaned.
tab high_vs_low wean_24hr
hist days_to_24hr_wean, by(high_vs_low)

gen censored_days_to_24hr_wean = days_to_24hr_wean

/* change comment for sensitivity analysis with noninformative censoring (discharge): */ 
replace censored_days_to_24hr_wean = 180 if missing(days_to_24hr_wean)
//replace sens_cens_days_to_24hr_wean = daystodischargefromrehab if missing(days_to_24hr_wean)
//note("Sensitivity analysis: Non-informative Censoring")

stset censored_days_to_24hr_wean, failure(wean_24hr==1)
sts test high_vs_low, logrank 
sts graph, by(high_vs_low) failure tmax(90) ///
 plotopts(lwidth(thick)) ///
 risktable(0(10)90, order(2 "Low Cervical" 1 "High Cervical") title("Number Eligible", size(medlarge)) size(medlarge)) ///
 xlabel(0(10)90, labsize(medlarge)) ///
 ylabel(,labsize(medlarge)) ///
 xtitle("Day of Rehab Stay", size(medlarge)) ///
 ytitle("Portion Weaned 24-Hr", size(medlarge)) ///
 text(0.9 5 "Log-Rank {it:p} = 0.053", placement(e) size(medlarge)) ///
 legend(order(1 "High Cervical (C4 & Above)" 2 "Low Cervical (C5 & Below)") position(4) ring(0) rows(2) size(medsmall)) ///
 title("Time 24-Hr Wean Occured", size(large))
graph export "Results and Figures/$S_DATE/Unadjusted KM 24-wean Level.png", as(png) name("Graph") replace
graph save "KM_24_wean.gph", replace
stcox high_vs_low

restore

/* "Survival" Until Decannulation */ 
preserve 
drop if partial_wean_at_admit == 2 //exclude those already partially weaned.
tab high_vs_low decannulate
hist daysfromadmissiontorehabtodecanu, by(high_vs_low)

gen censored_days_to_decan = daysfromadmissiontorehabtodecanu
replace censored_days_to_decan = 180 if missing(daysfromadmissiontorehabtodecanu)

stset censored_days_to_decan, failure(decannulate==1)
sts test high_vs_low, logrank 
sts graph, by(high_vs_low) failure tmax(90) ///
 plotopts(lwidth(thick)) ///
 text(0.7 5 "Log-Rank {it:p} = 0.009", placement(e) size(medlarge)) ///
 risktable(0(10)90, order(2 "Low Cervical" 1 "High Cervical") size(medlarge) title("Number Eligible", size(medlarge))) ///
 xlabel(0(10)90, labsize(medlarge)) ///
 ylabel(,labsize(medlarge)) /// 
 xtitle("Day of Rehab Stay", size(medlarge)) ///
 ytitle("Portion Decannulated", size(medlarge)) ///
 legend(order(1 "High Cervical (C4 & Above)" 2 "Low Cervical (C5 & Below)") position(10) ring(0) rows(2) size(medsmall)) ///
 text(0.03 30 "*One patient with a high CSCI decannulated at day 179", placement(e) size(small)) ///
 title("When Decannulation Occured", size(large)) 
graph export "Results and Figures/$S_DATE/Unadjusted KM Decannulation.png", as(png) name("Graph") replace
graph save "KM_decannulation.gph", replace
stcox high_vs_low
restore

//TODO: maybe add boxes around this?
graph combine KM_day_wean.gph KM_24_wean.gph KM_decannulation.gph, ///
	cols(1) /// 
	xsize(5) ysize(10)
graph export "Results and Figures/$S_DATE/Figure 2 - KMs for milestones.png", name("Graph") replace



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


/* Supplemental Table - does weaning status influence discharge location? */ 

ologit discharge_to ib1.weaning_outcome c.age, or
estimates store age_adj_wean_to_discharge
coefplot age_adj_wean_to_discharge, eform xscale(log) xline(1) xlabel(0.061 0.125 0.25 0.5 1 2 4 8 16 32 64 128) xscale(extend) xtitle("Odds Ratio of a higher level of care (LOC) discharge location" , size(small)) yscale(extend) ciopts(recast(rcap) lwidth(thick)) mlabel(string(@b,"%9.2f") + " [ " + string(@ll,"%9.2f") + " - " + string(@ul,"%9.2f") + " ] " + cond(@pval<.001, "***", cond(@pval<.01, "**", cond(@pval<.05, "*", "")))) mlabsize(medsmall) mlabposition(12) mlabgap(*1) scheme(white_tableau) text(5 0.11 "Favors Higher" "LOC Discharge" 5 9.9 "Favors Lower" "LOC Discharge", size(small) color(gs9)) headings(age = "{bf:Age} (per add'n year)" 2.weaning_outcome = "{bf:Weaning Outcome?} (vs 24h Vent)")
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

//TODO: change axis to years
stset time_to_censor_death, failure(death==1)
sts test discharge_to, logrank 

sts graph, by(discharge_to) tmax(2190) ///
 plotopts(lwidth(thick)) ///
 risktable(0(365)2190, order(1 "LTAC" 2 "SNF" 3 "Home with HH" 4 "Home") size(medlarge) title("Number Eligible", size(medlarge))) ///
 text(0.2 30 "Log-Rank {it:p} = 0.11", placement(e) size(medlarge)) ///
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
 text(0.2 30 "Log-Rank {it:p} = 0.11", placement(e) size(medlarge)) ///
 xlabel(0(365)2190, labsize(medlarge)) ///
 ylabel(, labsize(medlarge)) ///
 xtitle("Day From Rehab Admission", size(medlarge)) ///
 ytitle("Survival", size(medlarge)) ///
 legend(order(1 "Full Vent" 2 "Noct Vent" 3 "No Vent" 4 "Decannulated") position(4) ring(0) rows(2) size(medsmall)) ///
 title("Mortality by Discharge Ventilator Support", size(large))
graph export "Results and Figures/$S_DATE/Supp - KM Death by Wean.png", as(png) name("Graph") replace
graph save "KM_death_by_wean.gph", replace
stcox ib1.weaning_outcome

stcox ib1.weaning_outcome ib1.discharge_to 


graph combine KM_death_by_wean.gph KM_death_by_dispo.gph, ///
	cols(2) /// 
	xsize(10) ysize(5)
graph export "Results and Figures/$S_DATE/Figure 4 - KMs for death.png", name("Graph") replace



log close


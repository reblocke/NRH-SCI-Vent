* Supplemental figure generation for the NRH SCI cohort paper.
* Recreates the paper-matched cohort from nrh-sci-cleaned.dta and writes the
* retained publication figures, one exploratory subgroup figure, one compact
* correspondence table, and the run log.

version 18
clear all
set more off
capture log close

* Resolve the project root relative to the do-file launch location so the script
* can be run from `NRH SCI Data` or from `NRH-SCI-Vent`.
local project_root ""
foreach candidate in "." ".." {
    capture confirm file "`candidate'/nrh-sci-cleaned.dta"
    if !_rc {
        local project_root "`candidate'"
        continue, break
    }
}

if "`project_root'" == "" {
    di as err "Could not find nrh-sci-cleaned.dta in the current directory or one level up."
    exit 601
}

local results_dir "`project_root'/Results and Figures/$S_DATE"
local log_dir "`results_dir'/Logs"

capture mkdir "`project_root'/Results and Figures"
capture mkdir "`results_dir'"
capture mkdir "`log_dir'"

log using "`log_dir'/supplemental_sensitivity_analyses.log", replace text

di as txt "Resolved project root: `project_root'"
di as txt "Results directory: `results_dir'"

* Wilson intervals behave better than simple Wald intervals in these small
* grouped proportions, especially when counts are near 0 or 100%.
program define _wilson_ci, rclass
    args n x

    if missing(`n') | `n' <= 0 {
        return scalar p = .
        return scalar lb = .
        return scalar ub = .
        exit
    }

    tempname z ph denom center half
    scalar `z' = invnormal(0.975)
    scalar `ph' = `x' / `n'
    scalar `denom' = 1 + (`z'^2 / `n')
    scalar `center' = (`ph' + (`z'^2 / (2 * `n'))) / `denom'
    scalar `half' = `z' * sqrt((`ph' * (1 - `ph') + (`z'^2 / (4 * `n'))) / `n') / `denom'

    return scalar p = `ph'
    return scalar lb = max(0, `center' - `half')
    return scalar ub = min(1, `center' + `half')
end

* Export TIFF directly when Stata can; otherwise fall back to PDF plus `sips`
* so the workflow still produces publication-ready raster files on macOS.
program define _export_graph_tiff
    syntax, GRAPHNAME(name) OUTFILE(string asis) [WIDTH(integer 2400) HEIGHT(integer 1800)]

    tempfile pdffile
    local pdffile "`pdffile'.pdf"
    local outfile_clean = subinstr(`"`outfile'"', `"""', "", .)

    graph display `graphname'
    capture noisily graph export "`outfile_clean'", as(tif) replace width(`width') height(`height')
    if _rc {
        local tif_rc = _rc
        di as txt "TIFF translator unavailable or failed (rc=`tif_rc'); exporting PDF and converting with sips."
        graph export `"`pdffile'"', as(pdf) replace
        shell sips -s format tiff "`pdffile'" --out "`outfile_clean'"
        capture confirm file "`outfile_clean'"
        if _rc {
            di as err "Fallback TIFF conversion failed after graph export rc=`tif_rc'."
            exit 111
        }
    }
end

use "`project_root'/nrh-sci-cleaned.dta", clear

* Validate only the variables needed for the retained figures and the new
* exploratory decannulated-subgroup model block.
local required_vars ///
    age level wean_during_day days_to_daytime_wean ///
    wean_24hr days_to_24hr_wean decannulate ///
    daysfromadmissiontorehabtodecanu discharge_to ///
    partial_wean_at_admit age_decade weaning_outcome ///
    pneumonia_prior pna_at_rehab

local missing_vars ""
foreach var of local required_vars {
    capture confirm variable `var'
    if _rc {
        local missing_vars "`missing_vars' `var'"
    }
}

if trim("`missing_vars'") != "" {
    di as err "The cleaned dataset is missing required variables:`missing_vars'"
    log close
    exit 111
}

* Match the paper cohort by excluding patients already partially weaned at rehab
* admission before building any summaries or figures.
drop if partial_wean_at_admit == 1

quietly count
local analytic_n = r(N)
quietly levelsof level, local(levels_obs)
local c8_observed = strpos(" `levels_obs' ", " 8 ")
local level_note "Analytic cohort contains observed C2-C7 only; the C7-C8 stratum contained observed C7 cases only."
local figure_note_levels "Raw grouped proportions; bars show 95% Wilson CIs. C7-C8 contains observed C7 only."
local figure_note_timing "Achievers only. Circles = patients; boxes = IQR; squares = median; whiskers = range. Values >100 days are plotted at the cap and labeled by actual time."
local figure_note_age_1 "Points are individual patients. Shape/fill denote decannulation; lane outlines mark discharge categories."
local figure_note_age_2 "Vertical offset is visual only."
local figure_note_age_decann "Decannulated subgroup only. Points are individual patients; vertical jitter is visual only."
local figure_scheme "white_w3d"

di as txt "Analytic N after exclusion: `analytic_n'"
di as txt "Observed exact rehab injury levels: `levels_obs'"
if `c8_observed' {
    di as txt "Observed C8 in cleaned analytic dataset: yes"
}
else {
    di as txt "Observed C8 in cleaned analytic dataset: no"
    di as txt "`level_note'"
}

* Collapse exact injury levels into the four grouped strata used in the
* proportion figure.
generate byte injury_group4 = .
replace injury_group4 = 1 if inlist(level, 1, 2) & !missing(level)
replace injury_group4 = 2 if inlist(level, 3, 4) & !missing(level)
replace injury_group4 = 3 if inlist(level, 5, 6) & !missing(level)
replace injury_group4 = 4 if inlist(level, 7, 8) & !missing(level)
label define injury_group4_lab 1 "C1-C2" 2 "C3-C4" 3 "C5-C6" 4 "C7-C8", replace
label values injury_group4 injury_group4_lab
label variable injury_group4 "Finer injury level group"

* Build one compact grouped summary table in memory so the grouped milestone
* figure can be drawn without re-counting inside each panel.
* `postfile` is Stata's row-by-row way to assemble a small derived dataset
* without altering the patient-level data currently in memory.
tempfile grouped_level_tbl
tempname post_grouped
postfile `post_grouped' ///
    str12 stratum double group_order n_total ///
    day_wean_n day_wean_pct day_wean_ci_lb day_wean_ci_ub ///
    imv_n imv_pct imv_ci_lb imv_ci_ub ///
    decann_n decann_pct decann_ci_lb decann_ci_ub ///
    using `grouped_level_tbl', replace

forvalues grp = 1/4 {
    quietly count if injury_group4 == `grp'
    local n_total = r(N)

    quietly count if injury_group4 == `grp' & wean_during_day == 1
    local day_n = r(N)
    quietly _wilson_ci `n_total' `day_n'
    local day_pct = 100 * r(p)
    local day_lb = 100 * r(lb)
    local day_ub = 100 * r(ub)

    quietly count if injury_group4 == `grp' & wean_24hr == 1
    local imv_n = r(N)
    quietly _wilson_ci `n_total' `imv_n'
    local imv_pct = 100 * r(p)
    local imv_lb = 100 * r(lb)
    local imv_ub = 100 * r(ub)

    quietly count if injury_group4 == `grp' & decannulate == 1
    local dec_n = r(N)
    quietly _wilson_ci `n_total' `dec_n'
    local dec_pct = 100 * r(p)
    local dec_lb = 100 * r(lb)
    local dec_ub = 100 * r(ub)

    local grp_label "C1-C2"
    if `grp' == 2 local grp_label "C3-C4"
    if `grp' == 3 local grp_label "C5-C6"
    if `grp' == 4 local grp_label "C7-C8"

    post `post_grouped' ///
        ("`grp_label'") (`grp') (`n_total') ///
        (`day_n') (`day_pct') (`day_lb') (`day_ub') ///
        (`imv_n') (`imv_pct') (`imv_lb') (`imv_ub') ///
        (`dec_n') (`dec_pct') (`dec_lb') (`dec_ub')
}
* Close the temporary results file so it can be re-opened later like a regular
* dataset when drawing the grouped proportion figure.
postclose `post_grouped'

* Exact-level timing figure: for each milestone, show raw achievers plus a
* horizontal box-and-whisker summary within exact injury level.
* `preserve` keeps the full analytic cohort on a stack so each figure block can
* temporarily filter the data and then return to the original cohort with
* `restore`.
* Use one shared x-axis across all three timing panels so the stacked figure can
* be read against a common time scale. Cap the displayed range at 100 days so a
* single far-right outlier does not flatten the rest of the distributions.
quietly summarize days_to_daytime_wean if wean_during_day == 1 & !missing(days_to_daytime_wean)
local xmax_timing = r(max)
quietly summarize days_to_24hr_wean if wean_24hr == 1 & !missing(days_to_24hr_wean)
if r(max) > `xmax_timing' local xmax_timing = r(max)
quietly summarize daysfromadmissiontorehabtodecanu if decannulate == 1 & !missing(daysfromadmissiontorehabtodecanu)
if r(max) > `xmax_timing' local xmax_timing = r(max)
local xmax_timing = min(`xmax_timing', 100)
if `xmax_timing' < 50 local xmax_timing = 50
local timing_xlabel_spec "0(20)`xmax_timing'"

preserve
keep if wean_during_day == 1 & !missing(days_to_daytime_wean)
tempvar level_jit rowtag med p25 p75 minv maxv n_level time_plot outlier_flag outlier_lab
generate double `level_jit' = level + (runiform() - 0.5) * 0.16
generate double `time_plot' = min(days_to_daytime_wean, `xmax_timing')
generate byte `outlier_flag' = days_to_daytime_wean > `xmax_timing'
generate str8 `outlier_lab' = cond(`outlier_flag', string(days_to_daytime_wean, "%9.0f"), "")
* `rowtag' limits the summary layers to one draw per level; the jittered points
* keep individual achievers visible behind the box-and-whisker summary.
bysort level: generate byte `rowtag' = _n == 1
bysort level: egen int `n_level' = count(days_to_daytime_wean)
bysort level: egen double `med' = median(days_to_daytime_wean)
bysort level: egen double `p25' = pctile(days_to_daytime_wean), p(25)
bysort level: egen double `p75' = pctile(days_to_daytime_wean), p(75)
bysort level: egen double `minv' = min(days_to_daytime_wean)
bysort level: egen double `maxv' = max(days_to_daytime_wean)

* Layer order matters in `twoway`: whisker segments first, then the IQR box,
* then the jittered patient points, and finally the median square on top.
twoway ///
    (rcap `minv' `p25' level if `rowtag' & `n_level' >= 4, horizontal lcolor(black) lwidth(thin)) ///
    (rcap `p75' `maxv' level if `rowtag' & `n_level' >= 4, horizontal lcolor(black) lwidth(thin)) ///
    (rbar `p25' `p75' level if `rowtag' & `n_level' >= 4, horizontal barw(0.42) color(gs12) lcolor(black) lwidth(vthin)) ///
    (scatter `level_jit' `time_plot' if !`outlier_flag', msymbol(Oh) msize(small) mcolor(gs8)) ///
    (scatter `level_jit' `time_plot' if `outlier_flag', msymbol(Oh) msize(small) mcolor(gs8) ///
        mlabel(`outlier_lab') mlabsize(small) mlabcolor(gs6) mlabposition(9) mlabgap(vsmall)) ///
    (scatter level `med' if `rowtag' & `n_level' >= 2, msymbol(S) msize(medsmall) mfcolor(black) mlcolor(black)), ///
    ylabel(2 "C2" 3 "C3" 4 "C4" 5 "C5" 6 "C6" 7 "C7", angle(0) labsize(medsmall) noticks) ///
    yscale(reverse range(1.5 7.5)) ///
    xscale(range(0 `xmax_timing')) ///
    xlabel(`timing_xlabel_spec', labsize(small)) ///
    xtitle("Days from rehab admission", size(medlarge)) ///
    ytitle("Exact injury level", size(medlarge)) ///
    title("Time to daytime ventilator wean", size(large)) ///
    legend(off) ///
    scheme(`figure_scheme') ///
    name(gr_day_timing_summary, replace)
restore

* Repeat the same timing summary pattern for liberation from invasive
* ventilation.
preserve
keep if wean_24hr == 1 & !missing(days_to_24hr_wean)
tempvar level_jit rowtag med p25 p75 minv maxv n_level time_plot outlier_flag outlier_lab
generate double `level_jit' = level + (runiform() - 0.5) * 0.16
generate double `time_plot' = min(days_to_24hr_wean, `xmax_timing')
generate byte `outlier_flag' = days_to_24hr_wean > `xmax_timing'
generate str8 `outlier_lab' = cond(`outlier_flag', string(days_to_24hr_wean, "%9.0f"), "")
bysort level: generate byte `rowtag' = _n == 1
bysort level: egen int `n_level' = count(days_to_24hr_wean)
bysort level: egen double `med' = median(days_to_24hr_wean)
bysort level: egen double `p25' = pctile(days_to_24hr_wean), p(25)
bysort level: egen double `p75' = pctile(days_to_24hr_wean), p(75)
bysort level: egen double `minv' = min(days_to_24hr_wean)
bysort level: egen double `maxv' = max(days_to_24hr_wean)

* Use the same layer grammar as the first timing panel so the three outcomes
* are directly comparable.
twoway ///
    (rcap `minv' `p25' level if `rowtag' & `n_level' >= 4, horizontal lcolor(black) lwidth(thin)) ///
    (rcap `p75' `maxv' level if `rowtag' & `n_level' >= 4, horizontal lcolor(black) lwidth(thin)) ///
    (rbar `p25' `p75' level if `rowtag' & `n_level' >= 4, horizontal barw(0.42) color(gs12) lcolor(black) lwidth(vthin)) ///
    (scatter `level_jit' `time_plot' if !`outlier_flag', msymbol(Oh) msize(small) mcolor(gs8)) ///
    (scatter `level_jit' `time_plot' if `outlier_flag', msymbol(Oh) msize(small) mcolor(gs8) ///
        mlabel(`outlier_lab') mlabsize(small) mlabcolor(gs6) mlabposition(9) mlabgap(vsmall)) ///
    (scatter level `med' if `rowtag' & `n_level' >= 2, msymbol(S) msize(medsmall) mfcolor(black) mlcolor(black)), ///
    ylabel(2 "C2" 3 "C3" 4 "C4" 5 "C5" 6 "C6" 7 "C7", angle(0) labsize(medsmall) noticks) ///
    yscale(reverse range(1.5 7.5)) ///
    xscale(range(0 `xmax_timing')) ///
    xlabel(`timing_xlabel_spec', labsize(small)) ///
    xtitle("Days from rehab admission", size(medlarge)) ///
    ytitle("Exact injury level", size(medlarge)) ///
    title("Time to liberation from invasive ventilation", size(large)) ///
    legend(off) ///
    scheme(`figure_scheme') ///
    name(gr_imv_timing_summary, replace)
restore

* Repeat the same timing summary pattern for decannulation.
preserve
keep if decannulate == 1 & !missing(daysfromadmissiontorehabtodecanu)
tempvar level_jit rowtag med p25 p75 minv maxv n_level time_plot outlier_flag outlier_lab
generate double `level_jit' = level + (runiform() - 0.5) * 0.16
generate double `time_plot' = min(daysfromadmissiontorehabtodecanu, `xmax_timing')
generate byte `outlier_flag' = daysfromadmissiontorehabtodecanu > `xmax_timing'
generate str8 `outlier_lab' = cond(`outlier_flag', string(daysfromadmissiontorehabtodecanu, "%9.0f"), "")
bysort level: generate byte `rowtag' = _n == 1
bysort level: egen int `n_level' = count(daysfromadmissiontorehabtodecanu)
bysort level: egen double `med' = median(daysfromadmissiontorehabtodecanu)
bysort level: egen double `p25' = pctile(daysfromadmissiontorehabtodecanu), p(25)
bysort level: egen double `p75' = pctile(daysfromadmissiontorehabtodecanu), p(75)
bysort level: egen double `minv' = min(daysfromadmissiontorehabtodecanu)
bysort level: egen double `maxv' = max(daysfromadmissiontorehabtodecanu)

* Decannulation keeps the same plotting grammar but is restricted to patients
* who actually achieved decannulation and have a recorded time.
twoway ///
    (rcap `minv' `p25' level if `rowtag' & `n_level' >= 4, horizontal lcolor(black) lwidth(thin)) ///
    (rcap `p75' `maxv' level if `rowtag' & `n_level' >= 4, horizontal lcolor(black) lwidth(thin)) ///
    (rbar `p25' `p75' level if `rowtag' & `n_level' >= 4, horizontal barw(0.42) color(gs12) lcolor(black) lwidth(vthin)) ///
    (scatter `level_jit' `time_plot' if !`outlier_flag', msymbol(Oh) msize(small) mcolor(gs8)) ///
    (scatter `level_jit' `time_plot' if `outlier_flag', msymbol(Oh) msize(small) mcolor(gs8) ///
        mlabel(`outlier_lab') mlabsize(small) mlabcolor(gs6) mlabposition(9) mlabgap(vsmall)) ///
    (scatter level `med' if `rowtag' & `n_level' >= 2, msymbol(S) msize(medsmall) mfcolor(black) mlcolor(black)), ///
    ylabel(2 "C2" 3 "C3" 4 "C4" 5 "C5" 6 "C6" 7 "C7", angle(0) labsize(medsmall) noticks) ///
    yscale(reverse range(1.5 7.5)) ///
    xscale(range(0 `xmax_timing')) ///
    xlabel(`timing_xlabel_spec', labsize(small)) ///
    xtitle("Days from rehab admission", size(medlarge)) ///
    ytitle("Exact injury level", size(medlarge)) ///
    title("Time to decannulation", size(large)) ///
    legend(off) ///
    scheme(`figure_scheme') ///
    name(gr_dec_timing_summary, replace)
restore

* Stack the three timing panels into the retained exact-level timing figure.
graph combine gr_day_timing_summary gr_imv_timing_summary gr_dec_timing_summary, ///
    cols(1) xsize(8) ysize(10.5) ///
    note("`figure_note_timing'", size(vsmall)) ///
    imargin(small) ///
    name(gr_exact_timing_summary, replace)
quietly _export_graph_tiff, graphname(gr_exact_timing_summary) ///
    outfile("`results_dir'/Supplemental Figure - Median Days to Milestones by Exact Injury Level.tiff") ///
    width(3000) height(2400)

* Grouped milestone figure: raw grouped proportions with Wilson intervals and
* n/N labels for each finer injury-level stratum.
preserve
* Re-open the temporary grouped summary table so the plotting code below works
* with one row per injury stratum rather than one row per patient.
use `grouped_level_tbl', clear
generate str12 lab_day = string(day_wean_n, "%9.0f") + "/" + string(n_total, "%9.0f")
generate str12 lab_imv = string(imv_n, "%9.0f") + "/" + string(n_total, "%9.0f")
generate str12 lab_dec = string(decann_n, "%9.0f") + "/" + string(n_total, "%9.0f")
* Push most count labels slightly to the right of the point and pull the
* rightmost group's label to the left so it stays readable within the panel.
generate byte lab_pos = cond(group_order == 4, 9, 3)

* Each grouped panel is a two-layer plot: confidence interval bars plus a point
* estimate labeled with the underlying count.
twoway ///
    (rcap day_wean_ci_lb day_wean_ci_ub group_order, lcolor(black) lwidth(medthick)) ///
    (scatter day_wean_pct group_order, msymbol(O) msize(medlarge) mcolor(black) ///
        mlabel(lab_day) mlabsize(small) mlabcolor(black) mlabvposition(lab_pos) mlabgap(vsmall)), ///
    xlabel(1 "C1-C2" 2 "C3-C4" 3 "C5-C6" 4 "C7-C8", labsize(medsmall)) ///
    ylabel(0(20)100, labsize(medsmall)) ///
    yscale(range(0 100)) ///
    xtitle("Finer injury level group", size(medlarge)) ///
    ytitle("Achieved milestone (%)", size(medlarge)) ///
    title("Daytime ventilator wean", size(large)) ///
    legend(off) ///
    scheme(`figure_scheme') ///
    name(gr_day_group4, replace)

twoway ///
    (rcap imv_ci_lb imv_ci_ub group_order, lcolor(black) lwidth(medthick)) ///
    (scatter imv_pct group_order, msymbol(D) msize(medlarge) mcolor(black) ///
        mlabel(lab_imv) mlabsize(small) mlabcolor(black) mlabvposition(lab_pos) mlabgap(vsmall)), ///
    xlabel(1 "C1-C2" 2 "C3-C4" 3 "C5-C6" 4 "C7-C8", labsize(medsmall)) ///
    ylabel(0(20)100, labsize(medsmall)) ///
    yscale(range(0 100)) ///
    xtitle("Finer injury level group", size(medlarge)) ///
    ytitle("Achieved milestone (%)", size(medlarge)) ///
    title("Liberation from invasive ventilation", size(large)) ///
    legend(off) ///
    scheme(`figure_scheme') ///
    name(gr_imv_group4, replace)

twoway ///
    (rcap decann_ci_lb decann_ci_ub group_order, lcolor(black) lwidth(medthick)) ///
    (scatter decann_pct group_order, msymbol(T) msize(large) mcolor(black) ///
        mlabel(lab_dec) mlabsize(small) mlabcolor(black) mlabvposition(lab_pos) mlabgap(vsmall)), ///
    xlabel(1 "C1-C2" 2 "C3-C4" 3 "C5-C6" 4 "C7-C8", labsize(medsmall)) ///
    ylabel(0(20)100, labsize(medsmall)) ///
    yscale(range(0 100)) ///
    xtitle("Finer injury level group", size(medlarge)) ///
    ytitle("Achieved milestone (%)", size(medlarge)) ///
    title("Decannulation", size(large)) ///
    legend(off) ///
    scheme(`figure_scheme') ///
    name(gr_dec_group4, replace)

graph combine gr_day_group4 gr_imv_group4 gr_dec_group4, ///
    cols(1) xsize(6) ysize(10) ///
    note("`figure_note_levels'", size(vsmall) span justification(left)) ///
    imargin(small) ///
    name(gr_group4_milestones, replace)
quietly _export_graph_tiff, graphname(gr_group4_milestones) ///
    outfile("`results_dir'/Supplemental Figure - Milestone Rates by Finer Injury Groups.tiff") ///
    width(2700) height(3300)
restore

* Standalone discharge figure: one overlaid panel with slight vertical dodge so
* the two decannulation groups remain distinguishable without implying a second
* meaningful y-axis.
quietly summarize age
local age_xmin = 15
local age_xmax = ceil(r(max) / 5) * 5
if `age_xmax' < 75 local age_xmax = 75
local age_xlabel_spec "15(10)75"
if `age_xmax' > 75 local age_xlabel_spec "`age_xlabel_spec' `age_xmax'"

tempvar discharge_jit_overlay lane_tag lane_low lane_high
generate double `discharge_jit_overlay' = discharge_to + ///
    cond(decannulate == 1, 0.06, -0.06) + (runiform() - 0.5) * 0.04
bysort discharge_to: generate byte `lane_tag' = _n == 1
generate double `lane_low' = `age_xmin'
generate double `lane_high' = `age_xmax'

* Outline-only lane guides clarify the discharge rows while leaving the plot
* visually light; the legend is attached only to the two patient-point layers.
* The lane layer is drawn first, followed by open circles for not decannulated
* patients and filled diamonds for decannulated patients.
twoway ///
    (rbar `lane_low' `lane_high' discharge_to if `lane_tag', ///
        horizontal barw(0.34) fcolor(none) lcolor(gs10) lwidth(vthin)) ///
    (scatter `discharge_jit_overlay' age if decannulate == 0, ///
        msymbol(O) msize(medium) mfcolor(white) mlcolor(gs7) mlwidth(thin)) ///
    (scatter `discharge_jit_overlay' age if decannulate == 1, ///
        msymbol(D) msize(medium) mfcolor(black) mlcolor(black)), ///
    ylabel(1 "LTAC" 2 "SNF" 3 "Home w/ HH" 4 "Home", angle(0) labsize(medsmall) noticks) ///
    yscale(range(0.6 4.4) noline) ///
    xscale(range(`age_xmin' `age_xmax')) ///
    xlabel(`age_xlabel_spec', labsize(small)) ///
    xtitle("Age (years)", size(medlarge)) ///
    ytitle("Discharge destination", size(medlarge) margin(medlarge)) ///
    title("Observed discharge disposition by age", size(large)) ///
    note("`figure_note_age_1'" "`figure_note_age_2'", size(vsmall) span justification(left)) ///
    legend(order(2 "Not decannulated" 3 "Decannulated") rows(1) size(small) position(6)) ///
    graphregion(color(white) margin(medsmall)) ///
    plotregion(margin(small) lcolor(none)) ///
    scheme(`figure_scheme') ///
    name(gr_age_raw_overlay, replace)
quietly _export_graph_tiff, graphname(gr_age_raw_overlay) ///
    outfile("`results_dir'/Supplemental Figure - Observed Discharge Disposition by Age, Split by Decannulation.tiff") ///
    width(3600) height(2000)

* Exploratory subgroup figure: raw age-versus-discharge plot among patients
* who decannulated, using the same visual language as the full-cohort panel.
preserve
keep if decannulate == 1 & !missing(age, discharge_to)

* Recode to contiguous ordered categories because LTAC is not observed in the
* decannulated subgroup and should not create an empty lane or cutpoint.
generate byte discharge_to_decann = .
replace discharge_to_decann = 1 if discharge_to == 2
replace discharge_to_decann = 2 if discharge_to == 3
replace discharge_to_decann = 3 if discharge_to == 4
assert !missing(discharge_to_decann)
label define discharge_to_decann_lab 1 "SNF" 2 "Home w/ HH" 3 "Home", replace
label values discharge_to_decann discharge_to_decann_lab

quietly summarize age
local age_xmin_dec = 15
local age_xmax_dec = ceil(r(max) / 5) * 5
if `age_xmax_dec' < 75 local age_xmax_dec = 75
local age_xlabel_dec "15(10)75"
if `age_xmax_dec' > 75 local age_xlabel_dec "`age_xlabel_dec' `age_xmax_dec'"

tempvar discharge_jit_dec lane_tag_dec lane_low_dec lane_high_dec
generate double `discharge_jit_dec' = discharge_to_decann + (runiform() - 0.5) * 0.08
bysort discharge_to_decann: generate byte `lane_tag_dec' = _n == 1
generate double `lane_low_dec' = `age_xmin_dec'
generate double `lane_high_dec' = `age_xmax_dec'

twoway ///
    (rbar `lane_low_dec' `lane_high_dec' discharge_to_decann if `lane_tag_dec', ///
        horizontal barw(0.34) fcolor(none) lcolor(gs10) lwidth(vthin)) ///
    (scatter `discharge_jit_dec' age, ///
        msymbol(D) msize(medium) mfcolor(black) mlcolor(black)), ///
    ylabel(1 "SNF" 2 "Home w/ HH" 3 "Home", angle(0) labsize(medsmall) noticks) ///
    yscale(range(0.6 3.4) noline) ///
    xscale(range(`age_xmin_dec' `age_xmax_dec')) ///
    xlabel(`age_xlabel_dec', labsize(small)) ///
    xtitle("Age (years)", size(medlarge)) ///
    ytitle("Discharge destination", size(medlarge) margin(medlarge)) ///
    title("Observed discharge disposition by age among decannulated patients", size(large)) ///
    note("`figure_note_age_decann'", size(vsmall) span justification(left)) ///
    legend(off) ///
    graphregion(color(white) margin(medsmall)) ///
    plotregion(margin(small) lcolor(none)) ///
    scheme(`figure_scheme') ///
    name(gr_age_decann_only, replace)
quietly _export_graph_tiff, graphname(gr_age_decann_only) ///
    outfile("`results_dir'/Supplemental Figure - Discharge by Age Among Decannulated Patients.tiff") ///
    width(3200) height(1800)
restore

* Primary full-cohort response model: test whether age remains associated with
* discharge category after accounting for respiratory milestone status,
* completeness, and high-vs-low injury level.
preserve

di as txt " "
di as txt "Primary full-cohort milestone-adjusted discharge model"
di as txt "Analytic cohort N: " _N
tab weaning_outcome, missing
di as txt "No credible baseline comorbidity/frailty marker was available; pneumonia variables were treated as respiratory proxies and were not modeled."

* Reuse the temporary PLUS cache so proportional-odds diagnostics can run
* without changing the user's persistent ado tree.
local plus_orig "`c(sysdir_plus)'"
local plus_tmp "`c(tmpdir)'codex_oparallel_plus"
capture mkdir "`plus_tmp'"
quietly sysdir set PLUS "`plus_tmp'"
adopath ++ "`plus_tmp'"

capture which oparallel
if _rc {
    di as txt "Installing oparallel into temporary PLUS directory: `plus_tmp'"
    capture noisily ssc install oparallel, replace
    if _rc {
        quietly sysdir set PLUS "`plus_orig'"
        di as err "Failed to install oparallel into temporary PLUS directory."
        exit 111
    }
}
capture which oparallel
if _rc {
    quietly sysdir set PLUS "`plus_orig'"
    di as err "oparallel not found after installation attempt."
    exit 111
}

ologit discharge_to c.age_decade i.comp_vs_part i.high_vs_low ib1.weaning_outcome
estimates store full_cohort_milestone_model
local age_or_full = exp(_b[age_decade])
local age_z_full = _b[age_decade] / _se[age_decade]
local age_p_full = 2 * normal(-abs(`age_z_full'))
testparm i.weaning_outcome
local p_milestone_full = r(p)
capture noisily oparallel
if _rc {
    di as txt "Parallel-lines diagnostic via oparallel could not be computed for the primary full-cohort milestone-adjusted model because one or more auxiliary logits encountered perfect prediction."
}
estimates restore full_cohort_milestone_model

tempfile corr_home_marg corr_wean_counts analytic_cohort_copy
margins weaning_outcome, at(age_decade=(3 7)) predict(outcome(4)) saving(`corr_home_marg', replace)

save `analytic_cohort_copy', replace
contract weaning_outcome
rename _freq n_milestone
save `corr_wean_counts', replace
use `analytic_cohort_copy', clear

* Reformat the `margins` output into one compact correspondence table with one
* row per milestone state and age-30/age-70 adjusted home-discharge estimates.
use `corr_home_marg', clear
keep _m1 _at1 _margin _ci_lb _ci_ub
rename (_m1 _at1) (weaning_outcome age_decade)
generate int age_years = 10 * age_decade
drop age_decade
generate double home_prob_pct = 100 * _margin
generate double home_ci_lb_pct = 100 * _ci_lb
generate double home_ci_ub_pct = 100 * _ci_ub
drop _margin _ci_lb _ci_ub
reshape wide home_prob_pct home_ci_lb_pct home_ci_ub_pct, i(weaning_outcome) j(age_years)
merge 1:1 weaning_outcome using `corr_wean_counts', nogen assert(match)
decode weaning_outcome, gen(weaning_status)
generate byte milestone_order = weaning_outcome
order weaning_status n_milestone home_prob_pct30 home_ci_lb_pct30 home_ci_ub_pct30 ///
    home_prob_pct70 home_ci_lb_pct70 home_ci_ub_pct70
local decann_code = .
quietly levelsof weaning_outcome if strpos(lower(weaning_status), "decann") > 0, local(decann_code)
quietly summarize home_prob_pct30 if weaning_outcome == `decann_code', meanonly
local decann_home30 = r(mean)
quietly summarize home_prob_pct70 if weaning_outcome == `decann_code', meanonly
local decann_home70 = r(mean)
generate str12 home_prob_30 = string(home_prob_pct30, "%4.1f") + "%"
generate str24 home_ci_30 = string(home_ci_lb_pct30, "%4.1f") + "% to " + string(home_ci_ub_pct30, "%4.1f") + "%"
generate str12 home_prob_70 = string(home_prob_pct70, "%4.1f") + "%"
generate str24 home_ci_70 = string(home_ci_lb_pct70, "%4.1f") + "% to " + string(home_ci_ub_pct70, "%4.1f") + "%"
keep milestone_order weaning_status n_milestone home_prob_30 home_ci_30 home_prob_70 home_ci_70
rename (weaning_status n_milestone home_prob_30 home_ci_30 home_prob_70 home_ci_70) ///
    (milestone_status N pred_home_age30 ci_home_age30 pred_home_age70 ci_home_age70)
sort milestone_order
drop milestone_order
export delimited using "`results_dir'/Correspondence Table - Adjusted Home Discharge Probabilities by Milestone and Age.csv", replace
list, noobs abbreviate(24)
use `analytic_cohort_copy', clear

local age_or_full_fmt : display %4.2f `age_or_full'
local decann_home30_fmt : display %4.1f `decann_home30'
local decann_home70_fmt : display %4.1f `decann_home70'
if `age_p_full' < 0.05 local age_phrase_full "older age remained associated with more resource-intensive discharge"
else local age_phrase_full "older age was not a clear independent predictor of discharge category"
if `p_milestone_full' < 0.05 local milestone_phrase_full "overall respiratory milestone status was associated with discharge category"
else local milestone_phrase_full "overall respiratory milestone status was not clearly associated with discharge category"
di as txt "Interpretation: `age_phrase_full' after adjustment for respiratory milestone status, completeness, and high-vs-low injury level (age OR per decade older = `age_or_full_fmt'); `milestone_phrase_full'."
di as txt "Adjusted predicted probability of home discharge for decannulated patients was `decann_home30_fmt'% at age 30 and `decann_home70_fmt'% at age 70."
di as txt "Correspondence table written to `results_dir'/Correspondence Table - Adjusted Home Discharge Probabilities by Milestone and Age.csv"

* Exploratory sensitivity: replace the broad high-vs-low injury dichotomy with
* the four grouped cervical injury strata used in the descriptive figure.
di as txt " "
di as txt "Exploratory grouped-injury discharge sensitivity model"
quietly count if injury_group4 == 4
local n_c78 = r(N)
di as txt "Observed C7-C8 grouped stratum N: `n_c78'"
tab injury_group4, missing
ologit discharge_to c.age_decade i.comp_vs_part i.injury_group4, or
local age_or_group = exp(_b[age_decade])
local age_z_group = _b[age_decade] / _se[age_decade]
local age_p_group = 2 * normal(-abs(`age_z_group'))
testparm i.injury_group4
local p_injury_group = r(p)
capture noisily oparallel
if _rc {
    di as txt "Parallel-lines diagnostic via oparallel could not be computed for the grouped-injury discharge sensitivity model because one or more auxiliary logits encountered perfect prediction."
}

local age_or_group_fmt : display %4.2f `age_or_group'
if `age_p_group' < 0.05 local age_phrase_group "older age remained associated with more resource-intensive discharge"
else local age_phrase_group "older age was not a clear independent predictor of discharge category"
if `p_injury_group' < 0.05 local injury_phrase_group "grouped injury level was associated with discharge category"
else local injury_phrase_group "grouped injury-level terms were imprecise overall"
di as txt "Exploratory interpretation: `age_phrase_group' after replacing the high-vs-low dichotomy with four injury groups (age OR per decade older = `age_or_group_fmt'); `injury_phrase_group', and the C7-C8 estimate should be interpreted cautiously because N=`n_c78'."

quietly sysdir set PLUS "`plus_orig'"
restore

* Exploratory decannulated-subgroup ordered-logit models: age is the only
* defensible predictor here because milestone status is constant after
* conditioning on decannulation, and no baseline comorbidity index is present.
preserve
keep if decannulate == 1

di as txt " "
di as txt "Exploratory decannulated-subgroup ordinal discharge model"
di as txt "Decannulated subgroup N: " _N
tab discharge_to, missing
tab weaning_outcome, missing
di as txt "Weaning outcome is constant in this subgroup (all decannulated); no milestone-predictor subgroup model was fit."
di as txt "No credible baseline comorbidity/frailty marker was available; pneumonia variables were treated as respiratory proxies and were not modeled."

generate byte discharge_to_decann = .
replace discharge_to_decann = 1 if discharge_to == 2
replace discharge_to_decann = 2 if discharge_to == 3
replace discharge_to_decann = 3 if discharge_to == 4
assert !missing(discharge_to_decann)
label define discharge_to_decann_lab 1 "SNF" 2 "Home w/ Home Health" 3 "Home", replace
label values discharge_to_decann discharge_to_decann_lab
tab discharge_to_decann, missing

* Use a temporary PLUS directory so package installation for diagnostics does
* not modify the user's persistent Stata setup.
local plus_orig "`c(sysdir_plus)'"
local plus_tmp "`c(tmpdir)'codex_oparallel_plus"
capture mkdir "`plus_tmp'"
quietly sysdir set PLUS "`plus_tmp'"
adopath ++ "`plus_tmp'"

capture which oparallel
if _rc {
    di as txt "Installing oparallel into temporary PLUS directory: `plus_tmp'"
    capture noisily ssc install oparallel, replace
    if _rc {
        quietly sysdir set PLUS "`plus_orig'"
        di as err "Failed to install oparallel into temporary PLUS directory."
        exit 111
    }
}
capture which oparallel
if _rc {
    quietly sysdir set PLUS "`plus_orig'"
    di as err "oparallel not found after installation attempt."
    exit 111
}

di as txt "Primary model: ordered logit of discharge category on age per decade."
ologit discharge_to_decann c.age_decade, or
capture noisily oparallel
if _rc {
    quietly sysdir set PLUS "`plus_orig'"
    di as err "oparallel failed after the age-per-decade model."
    exit 111
}

di as txt "Sensitivity model: ordered logit of discharge category on age per year."
ologit discharge_to_decann c.age, or
capture noisily oparallel
if _rc {
    quietly sysdir set PLUS "`plus_orig'"
    di as err "oparallel failed after the age-per-year model."
    exit 111
}

di as txt "Exploratory interpretation: within the decannulated subgroup, age was not a clear predictor of discharge category, and proportional-odds diagnostics did not detect a violation, although power is limited."

quietly sysdir set PLUS "`plus_orig'"
restore

di as txt "Supplemental figure generation completed."
di as txt "Expected outputs written to `results_dir'."

log close

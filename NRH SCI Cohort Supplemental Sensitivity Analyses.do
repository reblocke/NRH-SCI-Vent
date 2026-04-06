version 18
clear all
set more off
capture log close

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

program define _model_diag_counts, rclass
    tempname b V conv
    matrix `b' = e(b)
    matrix `V' = e(V)

    local coefnames : colnames `b'
    local idx = 0
    local k_slopes = 0
    local omitted = 0

    foreach cname of local coefnames {
        local ++idx
        if substr("`cname'", 1, 1) == "/" {
            continue
        }
        if "`cname'" == "_cons" {
            continue
        }

        local ++k_slopes
        scalar `conv' = el(`V', `idx', `idx')
        if missing(`conv') | `conv' <= 0 {
            local ++omitted
        }
    }

    tempname estconv
    capture scalar `estconv' = e(converged)
    if !_rc {
        return scalar converged = scalar(`estconv')
    }
    else {
        return scalar converged = .
    }

    return scalar k_slopes = `k_slopes'
    return scalar omitted = `omitted'
end

program define _post_estimation_terms
    syntax, POSTNAME(name) MODELNAME(string asis) FAMILY(string asis) ESTTYPE(string asis) N(real)

    tempname b V coef var se est ll ul p
    matrix `b' = e(b)
    matrix `V' = e(V)
    local modelname_clean = subinstr(`"`modelname'"', `"""', "", .)
    local family_clean = subinstr(`"`family'"', `"""', "", .)
    local esttype_clean = subinstr(`"`esttype'"', `"""', "", .)

    local coefnames : colnames `b'
    local idx = 0

    foreach cname of local coefnames {
        local ++idx
        if substr("`cname'", 1, 1) == "/" {
            continue
        }
        if "`cname'" == "_cons" {
            continue
        }

        scalar `coef' = el(`b', 1, `idx')
        scalar `var' = el(`V', `idx', `idx')

        local row_status "ok"

        if missing(`var') | `var' <= 0 {
            scalar `est' = .
            scalar `ll' = .
            scalar `ul' = .
            scalar `p' = .
            local row_status "omitted"
        }
        else {
            scalar `se' = sqrt(`var')
            scalar `est' = exp(`coef')
            scalar `ll' = exp(`coef' - invnormal(0.975) * `se')
            scalar `ul' = exp(`coef' + invnormal(0.975) * `se')
            scalar `p' = 2 * normal(-abs(`coef' / `se'))
        }

        post `postname' ("`modelname_clean'") ("`family_clean'") ("`cname'") ("`esttype_clean'") ///
            (`n') (`est') (`ll') (`ul') (`p') ("`row_status'")
    }
end

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

local required_vars ///
    age age_decade level high_vs_low comp_vs_part reclass_on_arrival ///
    wean_during_day days_to_daytime_wean wean_24hr days_to_24hr_wean ///
    decannulate daysfromadmissiontorehabtodecanu weaning_outcome ///
    discharge_to home_vs_facility death daystodischargefromrehab ///
    partial_wean_at_admit

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

drop if partial_wean_at_admit == 1

quietly count
local analytic_n = r(N)
quietly levelsof level, local(levels_obs)
local c8_observed = strpos(" `levels_obs' ", " 8 ")
local level_note "Analytic cohort contains observed C2-C7 only; the C7-C8 stratum contained observed C7 cases only."
local small_sample_note "Small-sample support diagnostics are heuristic and summarize how much information each model has relative to its complexity; they should be read as caution flags rather than proof of overfitting."
local figure_note_levels "Supplemental sensitivity analysis. Points show observed proportions with 95% confidence intervals. `level_note'"
local figure_note_age_raw "Supplemental sensitivity analysis. Marker jitter is visual only."
local figure_note_age_model "Supplemental sensitivity analysis. Curves are model-based and should be interpreted cautiously given small subgroup counts."

di as txt "Analytic N after exclusion: `analytic_n'"
di as txt "Observed exact rehab injury levels: `levels_obs'"
if `c8_observed' {
    di as txt "Observed C8 in cleaned analytic dataset: yes"
}
else {
    di as txt "Observed C8 in cleaned analytic dataset: no"
    di as txt "`level_note'"
}

generate byte home_any = .
replace home_any = 1 if inlist(discharge_to, 3, 4) & !missing(discharge_to)
replace home_any = 0 if inlist(discharge_to, 1, 2) & !missing(discharge_to)
label define home_any_lab 0 "Facility (SNF, LTAC)" 1 "Home (+/- HH)", replace
label values home_any home_any_lab
label variable home_any "Discharge home or home with home health"

generate byte injury_group4 = .
replace injury_group4 = 1 if inlist(level, 1, 2) & !missing(level)
replace injury_group4 = 2 if inlist(level, 3, 4) & !missing(level)
replace injury_group4 = 3 if inlist(level, 5, 6) & !missing(level)
replace injury_group4 = 4 if inlist(level, 7, 8) & !missing(level)
label define injury_group4_lab 1 "C1-C2" 2 "C3-C4" 3 "C5-C6" 4 "C7-C8", replace
label values injury_group4 injury_group4_lab
label variable injury_group4 "Finer injury level group"

tempfile exact_level_tbl grouped_level_tbl level_notes_tbl ///
    ordered_model_tbl ordered_diag_tbl ordered_notes_tbl ///
    age_model_tbl age_diag_tbl age_notes_tbl ///
    decann_desc_tbl decann_model_tbl decann_notes_tbl

tempname post_exact post_grouped post_level_notes ///
    post_ordered_model post_ordered_diag post_ordered_notes ///
    post_age_model post_age_diag post_age_notes ///
    post_decann_desc post_decann_model post_decann_notes

postfile `post_exact' ///
    str12 stratum double group_order n_total ///
    day_wean_n day_wean_pct day_wean_ci_lb day_wean_ci_ub ///
    imv_n imv_pct imv_ci_lb imv_ci_ub ///
    decann_n decann_pct decann_ci_lb decann_ci_ub ///
    home_any_n home_any_pct home_any_ci_lb home_any_ci_ub ///
    day_wean_med day_wean_p25 day_wean_p75 ///
    imv_med imv_p25 imv_p75 ///
    decann_med decann_p25 decann_p75 ///
    using `exact_level_tbl', replace

postfile `post_grouped' ///
    str12 stratum double group_order n_total ///
    day_wean_n day_wean_pct day_wean_ci_lb day_wean_ci_ub ///
    imv_n imv_pct imv_ci_lb imv_ci_ub ///
    decann_n decann_pct decann_ci_lb decann_ci_ub ///
    home_any_n home_any_pct home_any_ci_lb home_any_ci_ub ///
    day_wean_med day_wean_p25 day_wean_p75 ///
    imv_med imv_p25 imv_p75 ///
    decann_med decann_p25 decann_p75 ///
    using `grouped_level_tbl', replace

postfile `post_level_notes' str244 note using `level_notes_tbl', replace

foreach lvl of local levels_obs {
    quietly count if level == `lvl'
    local n_total = r(N)

    quietly count if level == `lvl' & wean_during_day == 1
    local day_n = r(N)
    quietly _wilson_ci `n_total' `day_n'
    local day_pct = 100 * r(p)
    local day_lb = 100 * r(lb)
    local day_ub = 100 * r(ub)

    quietly count if level == `lvl' & wean_24hr == 1
    local imv_n = r(N)
    quietly _wilson_ci `n_total' `imv_n'
    local imv_pct = 100 * r(p)
    local imv_lb = 100 * r(lb)
    local imv_ub = 100 * r(ub)

    quietly count if level == `lvl' & decannulate == 1
    local dec_n = r(N)
    quietly _wilson_ci `n_total' `dec_n'
    local dec_pct = 100 * r(p)
    local dec_lb = 100 * r(lb)
    local dec_ub = 100 * r(ub)

    quietly count if level == `lvl' & home_any == 1
    local home_n = r(N)
    quietly _wilson_ci `n_total' `home_n'
    local home_pct = 100 * r(p)
    local home_lb = 100 * r(lb)
    local home_ub = 100 * r(ub)

    local day_med = .
    local day_p25 = .
    local day_p75 = .
    if `day_n' > 0 {
        quietly summarize days_to_daytime_wean if level == `lvl' & wean_during_day == 1, detail
        local day_med = r(p50)
        local day_p25 = r(p25)
        local day_p75 = r(p75)
    }

    local imv_med = .
    local imv_p25 = .
    local imv_p75 = .
    if `imv_n' > 0 {
        quietly summarize days_to_24hr_wean if level == `lvl' & wean_24hr == 1, detail
        local imv_med = r(p50)
        local imv_p25 = r(p25)
        local imv_p75 = r(p75)
    }

    local dec_med = .
    local dec_p25 = .
    local dec_p75 = .
    if `dec_n' > 0 {
        quietly summarize daysfromadmissiontorehabtodecanu if level == `lvl' & decannulate == 1, detail
        local dec_med = r(p50)
        local dec_p25 = r(p25)
        local dec_p75 = r(p75)
    }

    post `post_exact' ///
        ("C`lvl'") (`lvl') (`n_total') ///
        (`day_n') (`day_pct') (`day_lb') (`day_ub') ///
        (`imv_n') (`imv_pct') (`imv_lb') (`imv_ub') ///
        (`dec_n') (`dec_pct') (`dec_lb') (`dec_ub') ///
        (`home_n') (`home_pct') (`home_lb') (`home_ub') ///
        (`day_med') (`day_p25') (`day_p75') ///
        (`imv_med') (`imv_p25') (`imv_p75') ///
        (`dec_med') (`dec_p25') (`dec_p75')
}
postclose `post_exact'

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

    quietly count if injury_group4 == `grp' & home_any == 1
    local home_n = r(N)
    quietly _wilson_ci `n_total' `home_n'
    local home_pct = 100 * r(p)
    local home_lb = 100 * r(lb)
    local home_ub = 100 * r(ub)

    local day_med = .
    local day_p25 = .
    local day_p75 = .
    if `day_n' > 0 {
        quietly summarize days_to_daytime_wean if injury_group4 == `grp' & wean_during_day == 1, detail
        local day_med = r(p50)
        local day_p25 = r(p25)
        local day_p75 = r(p75)
    }

    local imv_med = .
    local imv_p25 = .
    local imv_p75 = .
    if `imv_n' > 0 {
        quietly summarize days_to_24hr_wean if injury_group4 == `grp' & wean_24hr == 1, detail
        local imv_med = r(p50)
        local imv_p25 = r(p25)
        local imv_p75 = r(p75)
    }

    local dec_med = .
    local dec_p25 = .
    local dec_p75 = .
    if `dec_n' > 0 {
        quietly summarize daysfromadmissiontorehabtodecanu if injury_group4 == `grp' & decannulate == 1, detail
        local dec_med = r(p50)
        local dec_p25 = r(p25)
        local dec_p75 = r(p75)
    }

    local grp_label "C1-C2"
    if `grp' == 2 local grp_label "C3-C4"
    if `grp' == 3 local grp_label "C5-C6"
    if `grp' == 4 local grp_label "C7-C8"

    post `post_grouped' ///
        ("`grp_label'") (`grp') (`n_total') ///
        (`day_n') (`day_pct') (`day_lb') (`day_ub') ///
        (`imv_n') (`imv_pct') (`imv_lb') (`imv_ub') ///
        (`dec_n') (`dec_pct') (`dec_lb') (`dec_ub') ///
        (`home_n') (`home_pct') (`home_lb') (`home_ub') ///
        (`day_med') (`day_p25') (`day_p75') ///
        (`imv_med') (`imv_p25') (`imv_p75') ///
        (`dec_med') (`dec_p25') (`dec_p75')
}
postclose `post_grouped'

post `post_level_notes' ("Analytic N after paper-matched exclusion: `analytic_n'")
post `post_level_notes' ("Observed exact rehab injury levels: `levels_obs'")
if `c8_observed' {
    post `post_level_notes' ("C8 was observed in the cleaned analytic dataset.")
}
else {
    post `post_level_notes' ("C8 was not observed in the cleaned analytic dataset.")
    post `post_level_notes' ("`level_note'")
}
post `post_level_notes' ("These descriptive summaries are supplemental sensitivity analyses and should be interpreted alongside, not instead of, the main paper analyses.")
postclose `post_level_notes'

preserve
use `exact_level_tbl', clear
export excel using "`results_dir'/Supplemental Table - Milestone Rates by Injury Level.xlsx", ///
    sheet("Exact injury level") firstrow(variables) replace
restore

preserve
use `grouped_level_tbl', clear
export excel using "`results_dir'/Supplemental Table - Milestone Rates by Injury Level.xlsx", ///
    sheet("Grouped injury level") firstrow(variables) sheetreplace
restore

preserve
use `level_notes_tbl', clear
export excel using "`results_dir'/Supplemental Table - Milestone Rates by Injury Level.xlsx", ///
    sheet("Notes") firstrow(variables) sheetreplace
restore

preserve
use `grouped_level_tbl', clear
generate str12 lab_day = string(day_wean_n, "%9.0f") + "/" + string(n_total, "%9.0f")
generate str12 lab_imv = string(imv_n, "%9.0f") + "/" + string(n_total, "%9.0f")
generate str12 lab_dec = string(decann_n, "%9.0f") + "/" + string(n_total, "%9.0f")

twoway ///
    (rcap day_wean_ci_lb day_wean_ci_ub group_order, lcolor(black) lwidth(medthick)) ///
    (scatter day_wean_pct group_order, msymbol(O) msize(medlarge) mcolor(black) ///
        mlabel(lab_day) mlabsize(small) mlabcolor(black) mlabposition(12)), ///
    xlabel(1 "C1-C2" 2 "C3-C4" 3 "C5-C6" 4 "C7-C8", labsize(medsmall)) ///
    ylabel(0(20)100, labsize(medsmall)) ///
    yscale(range(0 100)) ///
    xtitle("Finer injury level group", size(medlarge)) ///
    ytitle("Achieved milestone (%)", size(medlarge)) ///
    title("Daytime ventilator wean", size(large)) ///
    legend(off) ///
    scheme(s1mono) ///
    name(gr_day_group4, replace)

twoway ///
    (rcap imv_ci_lb imv_ci_ub group_order, lcolor(black) lwidth(medthick)) ///
    (scatter imv_pct group_order, msymbol(D) msize(medlarge) mcolor(black) ///
        mlabel(lab_imv) mlabsize(small) mlabcolor(black) mlabposition(12)), ///
    xlabel(1 "C1-C2" 2 "C3-C4" 3 "C5-C6" 4 "C7-C8", labsize(medsmall)) ///
    ylabel(0(20)100, labsize(medsmall)) ///
    yscale(range(0 100)) ///
    xtitle("Finer injury level group", size(medlarge)) ///
    ytitle("Achieved milestone (%)", size(medlarge)) ///
    title("Liberation from invasive ventilation", size(large)) ///
    legend(off) ///
    scheme(s1mono) ///
    name(gr_imv_group4, replace)

twoway ///
    (rcap decann_ci_lb decann_ci_ub group_order, lcolor(black) lwidth(medthick)) ///
    (scatter decann_pct group_order, msymbol(T) msize(large) mcolor(black) ///
        mlabel(lab_dec) mlabsize(small) mlabcolor(black) mlabposition(12)), ///
    xlabel(1 "C1-C2" 2 "C3-C4" 3 "C5-C6" 4 "C7-C8", labsize(medsmall)) ///
    ylabel(0(20)100, labsize(medsmall)) ///
    yscale(range(0 100)) ///
    xtitle("Finer injury level group", size(medlarge)) ///
    ytitle("Achieved milestone (%)", size(medlarge)) ///
    title("Decannulated", size(large)) ///
    legend(off) ///
    scheme(s1mono) ///
    name(gr_dec_group4, replace)

graph combine gr_day_group4 gr_imv_group4 gr_dec_group4, ///
    cols(1) xsize(6) ysize(10) ///
    note("`figure_note_levels'", size(vsmall)) ///
    imargin(small) ///
    name(gr_group4_milestones, replace)

graph save gr_group4_milestones ///
    "`results_dir'/Supplemental Figure - Milestone Rates by Finer Injury Groups.gph", replace
quietly _export_graph_tiff, graphname(gr_group4_milestones) ///
    outfile("`results_dir'/Supplemental Figure - Milestone Rates by Finer Injury Groups.tiff") ///
    width(2700) height(3300)
restore

postfile `post_ordered_model' ///
    str80 model_name str12 family str80 term str8 estimate_type ///
    double n estimate ci_lb ci_ub p_value str120 status ///
    using `ordered_model_tbl', replace

postfile `post_ordered_diag' ///
    str80 model_name str12 family double n cat1 cat2 cat3 cat4 failures compete ///
    k_slopes omitted_slopes converged heuristic_ratio ///
    str12 overfit_flag str120 status ///
    using `ordered_diag_tbl', replace

postfile `post_ordered_notes' str244 note using `ordered_notes_tbl', replace

capture noisily ologit weaning_outcome c.level c.age_decade i.comp_vs_part, or
if _rc {
    local fail_msg "weaning_outcome ordered-level model omitted (rc=`_rc')"
    di as err "`fail_msg'"
    post `post_ordered_model' ("Weaning outcome ordered-level sensitivity") ("ologit") ///
        ("<model omitted>") ("OR") (`analytic_n') (.) (.) (.) (.) ("`fail_msg'")
    post `post_ordered_diag' ("Weaning outcome ordered-level sensitivity") ("ologit") ///
        (`analytic_n') (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) ("High concern") ("`fail_msg'")
}
else {
    quietly _post_estimation_terms, postname(`post_ordered_model') ///
        modelname("Weaning outcome ordered-level sensitivity") ///
        family("ologit") esttype("OR") n(`analytic_n')

    quietly _model_diag_counts
    local k_slopes = r(k_slopes)
    local omitted = r(omitted)
    local converged = r(converged)

    local min_cat = .
    forvalues cat = 1/4 {
        quietly count if e(sample) & weaning_outcome == `cat'
        local c`cat' = r(N)
        if missing(`min_cat') | `c`cat'' < `min_cat' {
            local min_cat = `c`cat''
        }
    }

    local heuristic_ratio = .
    if `k_slopes' > 0 {
        local heuristic_ratio = `min_cat' / `k_slopes'
    }

    local overfit_flag "Lower concern"
    if `heuristic_ratio' < 10 {
        local overfit_flag "Moderate concern"
    }
    if `heuristic_ratio' < 5 | `omitted' > 0 | `converged' != 1 {
        local overfit_flag "High concern"
    }

    post `post_ordered_diag' ("Weaning outcome ordered-level sensitivity") ("ologit") ///
        (`analytic_n') (`c1') (`c2') (`c3') (`c4') (.) (.) ///
        (`k_slopes') (`omitted') (`converged') (`heuristic_ratio') ///
        ("`overfit_flag'") ("ok")
    di as txt "Weaning outcome ordered-level sensitivity heuristic ratio: `heuristic_ratio' (`overfit_flag')."
}

capture noisily ologit discharge_to c.level c.age_decade i.comp_vs_part, or
if _rc {
    local fail_msg "discharge_to ordered-level model omitted (rc=`_rc')"
    di as err "`fail_msg'"
    post `post_ordered_model' ("Discharge ordered-level sensitivity") ("ologit") ///
        ("<model omitted>") ("OR") (`analytic_n') (.) (.) (.) (.) ("`fail_msg'")
    post `post_ordered_diag' ("Discharge ordered-level sensitivity") ("ologit") ///
        (`analytic_n') (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) ("High concern") ("`fail_msg'")
}
else {
    quietly _post_estimation_terms, postname(`post_ordered_model') ///
        modelname("Discharge ordered-level sensitivity") ///
        family("ologit") esttype("OR") n(`analytic_n')

    quietly _model_diag_counts
    local k_slopes = r(k_slopes)
    local omitted = r(omitted)
    local converged = r(converged)

    local min_cat = .
    forvalues cat = 1/4 {
        quietly count if e(sample) & discharge_to == `cat'
        local c`cat' = r(N)
        if missing(`min_cat') | `c`cat'' < `min_cat' {
            local min_cat = `c`cat''
        }
    }

    local heuristic_ratio = .
    if `k_slopes' > 0 {
        local heuristic_ratio = `min_cat' / `k_slopes'
    }

    local overfit_flag "Lower concern"
    if `heuristic_ratio' < 10 {
        local overfit_flag "Moderate concern"
    }
    if `heuristic_ratio' < 5 | `omitted' > 0 | `converged' != 1 {
        local overfit_flag "High concern"
    }

    post `post_ordered_diag' ("Discharge ordered-level sensitivity") ("ologit") ///
        (`analytic_n') (`c1') (`c2') (`c3') (`c4') (.) (.) ///
        (`k_slopes') (`omitted') (`converged') (`heuristic_ratio') ///
        ("`overfit_flag'") ("ok")
    di as txt "Discharge ordered-level sensitivity heuristic ratio: `heuristic_ratio' (`overfit_flag')."
}

local milestones "daytime imv decann"
foreach milestone of local milestones {
    preserve

    generate byte event_type = .
    generate double time_to_event = .

    local model_name ""
    if "`milestone'" == "daytime" {
        local model_name "Daytime wean Fine-Gray ordered-level sensitivity"
        replace event_type = 1 if wean_during_day == 1
        replace event_type = 2 if missing(event_type) & !missing(daystodischargefromrehab) & ///
            (missing(days_to_daytime_wean) | daystodischargefromrehab < days_to_daytime_wean)
        replace time_to_event = days_to_daytime_wean if wean_during_day == 1 & !missing(days_to_daytime_wean)
        replace time_to_event = daystodischargefromrehab if event_type == 2
    }
    else if "`milestone'" == "imv" {
        local model_name "Liberation from IMV Fine-Gray ordered-level sensitivity"
        replace event_type = 1 if wean_24hr == 1
        replace event_type = 2 if missing(event_type) & !missing(daystodischargefromrehab) & ///
            (missing(days_to_24hr_wean) | daystodischargefromrehab < days_to_24hr_wean)
        replace time_to_event = days_to_24hr_wean if wean_24hr == 1 & !missing(days_to_24hr_wean)
        replace time_to_event = daystodischargefromrehab if event_type == 2
    }
    else if "`milestone'" == "decann" {
        local model_name "Decannulation Fine-Gray ordered-level sensitivity"
        replace event_type = 1 if decannulate == 1
        replace event_type = 2 if missing(event_type) & !missing(daystodischargefromrehab) & ///
            (missing(daysfromadmissiontorehabtodecanu) | daystodischargefromrehab < daysfromadmissiontorehabtodecanu)
        replace time_to_event = daysfromadmissiontorehabtodecanu if decannulate == 1 & !missing(daysfromadmissiontorehabtodecanu)
        replace time_to_event = daystodischargefromrehab if event_type == 2
    }

    quietly count if !missing(time_to_event)
    local model_n = r(N)

    capture noisily stset time_to_event, failure(event_type == 1)
    if _rc {
        local fail_msg "`model_name' stset omitted (rc=`_rc')"
        di as err "`fail_msg'"
        post `post_ordered_model' ("`model_name'") ("stcrreg") ///
            ("<model omitted>") ("SHR") (`model_n') (.) (.) (.) (.) ("`fail_msg'")
        post `post_ordered_diag' ("`model_name'") ("stcrreg") ///
            (`model_n') (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) ("High concern") ("`fail_msg'")
        restore
        continue
    }

    capture noisily stcrreg c.level c.age_decade i.comp_vs_part, compete(event_type == 2)
    if _rc {
        local fail_msg "`model_name' omitted (rc=`_rc')"
        di as err "`fail_msg'"
        post `post_ordered_model' ("`model_name'") ("stcrreg") ///
            ("<model omitted>") ("SHR") (`model_n') (.) (.) (.) (.) ("`fail_msg'")
        quietly count if event_type == 1
        local failures = r(N)
        quietly count if event_type == 2
        local compete = r(N)
        post `post_ordered_diag' ("`model_name'") ("stcrreg") ///
            (`model_n') (.) (.) (.) (.) (`failures') (`compete') ///
            (.) (.) (.) (.) ("High concern") ("`fail_msg'")
        restore
        continue
    }

    quietly _post_estimation_terms, postname(`post_ordered_model') ///
        modelname("`model_name'") family("stcrreg") esttype("SHR") n(`model_n')

    quietly _model_diag_counts
    local k_slopes = r(k_slopes)
    local omitted = r(omitted)
    local converged = r(converged)

    quietly count if e(sample) & event_type == 1
    local failures = r(N)
    quietly count if e(sample) & event_type == 2
    local compete = r(N)

    local heuristic_ratio = .
    if `k_slopes' > 0 {
        local heuristic_ratio = `failures' / `k_slopes'
    }

    local overfit_flag "Lower concern"
    if `heuristic_ratio' < 10 {
        local overfit_flag "Moderate concern"
    }
    if `heuristic_ratio' < 5 | `omitted' > 0 | `converged' != 1 {
        local overfit_flag "High concern"
    }

    post `post_ordered_diag' ("`model_name'") ("stcrreg") ///
        (`model_n') (.) (.) (.) (.) (`failures') (`compete') ///
        (`k_slopes') (`omitted') (`converged') (`heuristic_ratio') ///
        ("`overfit_flag'") ("ok")
    di as txt "`model_name' heuristic ratio: `heuristic_ratio' (`overfit_flag')."

    restore
}

post `post_ordered_notes' ("Ordered-level models treat injury level as an ordered exploratory sensitivity term; they do not imply proven linearity.")
post `post_ordered_notes' ("Fine-Gray sensitivity models mirror the paper's competing-risk structure but replace i.high_vs_low with c.level.")
post `post_ordered_notes' ("`small_sample_note'")
post `post_ordered_notes' ("For ordered-logit models, support is summarized as the smallest outcome-category count divided by the number of free slope parameters; for Fine-Gray models, support is summarized as the number of target failures divided by the number of free slope parameters.")
postclose `post_ordered_model'
postclose `post_ordered_diag'
postclose `post_ordered_notes'

preserve
use `ordered_model_tbl', clear
keep if family == "ologit"
export excel using "`results_dir'/Supplemental Table - Ordered Level Sensitivity Models.xlsx", ///
    sheet("Ordered-logit models") firstrow(variables) replace
restore

preserve
use `ordered_model_tbl', clear
keep if family == "stcrreg"
export excel using "`results_dir'/Supplemental Table - Ordered Level Sensitivity Models.xlsx", ///
    sheet("Fine-Gray models") firstrow(variables) sheetreplace
restore

preserve
use `ordered_diag_tbl', clear
export excel using "`results_dir'/Supplemental Table - Ordered Level Sensitivity Models.xlsx", ///
    sheet("Diagnostics") firstrow(variables) sheetreplace
restore

preserve
use `ordered_notes_tbl', clear
export excel using "`results_dir'/Supplemental Table - Ordered Level Sensitivity Models.xlsx", ///
    sheet("Notes") firstrow(variables) sheetreplace
restore

postfile `post_age_model' ///
    str80 model_name str12 family str80 term str8 estimate_type ///
    double n estimate ci_lb ci_ub p_value str120 status ///
    using `age_model_tbl', replace

postfile `post_age_diag' ///
    str80 model_name str12 family double n events nonevents cat1 cat2 cat3 cat4 ///
    k_slopes omitted_slopes converged heuristic_ratio linktest_p ///
    str12 overfit_flag str120 status ///
    using `age_diag_tbl', replace

postfile `post_age_notes' str244 note using `age_notes_tbl', replace

capture noisily ologit discharge_to c.age_decade i.weaning_outcome i.high_vs_low i.comp_vs_part, or
if _rc {
    local fail_msg "Milestone-adjusted discharge ologit omitted (rc=`_rc')"
    di as err "`fail_msg'"
    post `post_age_model' ("Milestone-adjusted discharge ordinal model") ("ologit") ///
        ("<model omitted>") ("OR") (`analytic_n') (.) (.) (.) (.) ("`fail_msg'")
    post `post_age_diag' ("Milestone-adjusted discharge ordinal model") ("ologit") ///
        (`analytic_n') (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) ("High concern") ("`fail_msg'")
}
else {
    quietly _post_estimation_terms, postname(`post_age_model') ///
        modelname("Milestone-adjusted discharge ordinal model") ///
        family("ologit") esttype("OR") n(`analytic_n')

    quietly _model_diag_counts
    local k_slopes = r(k_slopes)
    local omitted = r(omitted)
    local converged = r(converged)

    local min_cat = .
    forvalues cat = 1/4 {
        quietly count if e(sample) & discharge_to == `cat'
        local c`cat' = r(N)
        if missing(`min_cat') | `c`cat'' < `min_cat' {
            local min_cat = `c`cat''
        }
    }

    local heuristic_ratio = .
    if `k_slopes' > 0 {
        local heuristic_ratio = `min_cat' / `k_slopes'
    }

    local overfit_flag "Lower concern"
    if `heuristic_ratio' < 10 {
        local overfit_flag "Moderate concern"
    }
    if `heuristic_ratio' < 5 | `omitted' > 0 | `converged' != 1 {
        local overfit_flag "High concern"
    }

    post `post_age_diag' ("Milestone-adjusted discharge ordinal model") ("ologit") ///
        (`analytic_n') (.) (.) (`c1') (`c2') (`c3') (`c4') ///
        (`k_slopes') (`omitted') (`converged') (`heuristic_ratio') (.) ///
        ("`overfit_flag'") ("ok")
    di as txt "Milestone-adjusted discharge ologit heuristic ratio: `heuristic_ratio' (`overfit_flag')."
}

local interaction_fit_ok 0
local interaction_status "ok"

capture noisily logit home_any c.age_decade##i.decannulate i.high_vs_low i.comp_vs_part, or
if _rc {
    local interaction_status "Age-by-decannulation interaction logit omitted (rc=`_rc')"
    di as err "`interaction_status'"
    post `post_age_model' ("Age x decannulation home-any model") ("logit") ///
        ("<model omitted>") ("OR") (`analytic_n') (.) (.) (.) (.) ("`interaction_status'")
    post `post_age_diag' ("Age x decannulation home-any model") ("logit") ///
        (`analytic_n') (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) ("High concern") ("`interaction_status'")
}
else {
    local interaction_fit_ok 1
    estimates store age_home_interaction_model
    quietly _post_estimation_terms, postname(`post_age_model') ///
        modelname("Age x decannulation home-any model") ///
        family("logit") esttype("OR") n(`analytic_n')

    quietly _model_diag_counts
    local k_slopes = r(k_slopes)
    local omitted = r(omitted)
    local converged = r(converged)

    quietly count if e(sample) & home_any == 1
    local events = r(N)
    quietly count if e(sample) & home_any == 0
    local nonevents = r(N)

    local heuristic_ratio = .
    if `k_slopes' > 0 {
        local heuristic_ratio = min(`events', `nonevents') / `k_slopes'
    }

    local link_p = .
    local link_note = "not run"
    tempname link_b link_se
    capture noisily linktest, nolog
    if _rc {
        local link_note = "linktest failed (rc=`_rc')"
        di as txt "`link_note'"
    }
    else {
        capture scalar `link_b' = _b[_hatsq]
        if !_rc {
            scalar `link_se' = _se[_hatsq]
            local link_p = 2 * normal(-abs(`link_b' / `link_se'))
            if `link_p' < 0.05 {
                local link_note = "possible misspecification"
            }
            else {
                local link_note = "no strong misspecification signal"
            }
        }
        else {
            local link_note = "linktest completed; _hatsq unavailable"
        }
        estimates restore age_home_interaction_model
    }

    local overfit_flag "Lower concern"
    if `heuristic_ratio' < 10 {
        local overfit_flag "Moderate concern"
    }
    if `heuristic_ratio' < 5 | `omitted' > 0 | `converged' != 1 {
        local overfit_flag "High concern"
    }

    post `post_age_diag' ("Age x decannulation home-any model") ("logit") ///
        (`analytic_n') (`events') (`nonevents') (.) (.) (.) (.) ///
        (`k_slopes') (`omitted') (`converged') (`heuristic_ratio') (`link_p') ///
        ("`overfit_flag'") ("`link_note'")
    di as txt "Age x decannulation interaction heuristic ratio: `heuristic_ratio' (`overfit_flag')."
}

postfile `post_decann_desc' ///
    str24 subgroup double n_total home_any_n facility_n ///
    age_med_home age_p25_home age_p75_home ///
    age_med_facility age_p25_facility age_p75_facility ///
    using `decann_desc_tbl', replace

postfile `post_decann_model' ///
    str80 model_name str12 family str80 term str8 estimate_type ///
    double n estimate ci_lb ci_ub p_value str120 status ///
    using `decann_model_tbl', replace

postfile `post_decann_notes' str244 note using `decann_notes_tbl', replace

quietly count if decannulate == 1
local decann_n = r(N)
quietly count if decannulate == 1 & home_any == 1
local decann_home_n = r(N)
quietly count if decannulate == 1 & home_any == 0
local decann_fac_n = r(N)

local age_med_home = .
local age_p25_home = .
local age_p75_home = .
if `decann_home_n' > 0 {
    quietly summarize age if decannulate == 1 & home_any == 1, detail
    local age_med_home = r(p50)
    local age_p25_home = r(p25)
    local age_p75_home = r(p75)
}

local age_med_fac = .
local age_p25_fac = .
local age_p75_fac = .
if `decann_fac_n' > 0 {
    quietly summarize age if decannulate == 1 & home_any == 0, detail
    local age_med_fac = r(p50)
    local age_p25_fac = r(p25)
    local age_p75_fac = r(p75)
}

post `post_decann_desc' ("All decannulated") (`decann_n') (`decann_home_n') (`decann_fac_n') ///
    (`age_med_home') (`age_p25_home') (`age_p75_home') ///
    (`age_med_fac') (`age_p25_fac') (`age_p75_fac')
postclose `post_decann_desc'

quietly count if decannulate == 1 & home_any == 1 & high_vs_low == 1
local hh_home = r(N)
quietly count if decannulate == 1 & home_any == 0 & high_vs_low == 1
local hh_fac = r(N)
quietly count if decannulate == 1 & home_any == 1 & high_vs_low == 2
local ll_home = r(N)
quietly count if decannulate == 1 & home_any == 0 & high_vs_low == 2
local ll_fac = r(N)
quietly count if decannulate == 1 & home_any == 1 & comp_vs_part == 1
local comp_home = r(N)
quietly count if decannulate == 1 & home_any == 0 & comp_vs_part == 1
local comp_fac = r(N)
quietly count if decannulate == 1 & home_any == 1 & comp_vs_part == 2
local part_home = r(N)
quietly count if decannulate == 1 & home_any == 0 & comp_vs_part == 2
local part_fac = r(N)

local subgroup_model_ok 1
if `decann_home_n' < 5 | `decann_fac_n' < 5 local subgroup_model_ok 0
if `hh_home' == 0 | `hh_fac' == 0 | `ll_home' == 0 | `ll_fac' == 0 local subgroup_model_ok 0
if `comp_home' == 0 | `comp_fac' == 0 | `part_home' == 0 | `part_fac' == 0 local subgroup_model_ok 0

if !`subgroup_model_ok' {
    local subgroup_msg "Model omitted because subgroup cell counts were too sparse, raising concern for quasi-complete separation."
    di as txt "`subgroup_msg'"
    post `post_decann_model' ("Decannulated-only home-any model") ("logit") ///
        ("<model omitted>") ("OR") (`decann_n') (.) (.) (.) (.) ("`subgroup_msg'")
    post `post_age_diag' ("Decannulated-only home-any model") ("logit") ///
        (`decann_n') (`decann_home_n') (`decann_fac_n') (.) (.) (.) (.) ///
        (.) (.) (.) (.) (.) ("High concern") ("`subgroup_msg'")
    post `post_decann_notes' ("Interpret this subgroup descriptively only because the available information is limited.")
    post `post_decann_notes' ("`subgroup_msg'")
}
else {
    capture noisily logit home_any c.age_decade i.high_vs_low i.comp_vs_part if decannulate == 1, or
    if _rc {
        local subgroup_msg "Decannulated-only model omitted (rc=`_rc')"
        di as err "`subgroup_msg'"
        post `post_decann_model' ("Decannulated-only home-any model") ("logit") ///
            ("<model omitted>") ("OR") (`decann_n') (.) (.) (.) (.) ("`subgroup_msg'")
        post `post_age_diag' ("Decannulated-only home-any model") ("logit") ///
            (`decann_n') (`decann_home_n') (`decann_fac_n') (.) (.) (.) (.) ///
            (.) (.) (.) (.) (.) ("High concern") ("`subgroup_msg'")
        post `post_decann_notes' ("Interpret this subgroup descriptively only because the fitted model was unstable.")
    }
    else {
        quietly _post_estimation_terms, postname(`post_decann_model') ///
            modelname("Decannulated-only home-any model") ///
            family("logit") esttype("OR") n(`decann_n')

        quietly _model_diag_counts
        local k_slopes = r(k_slopes)
        local omitted = r(omitted)
        local converged = r(converged)
        local heuristic_ratio = .
        if `k_slopes' > 0 {
            local heuristic_ratio = min(`decann_home_n', `decann_fac_n') / `k_slopes'
        }

        local overfit_flag "Lower concern"
        if `heuristic_ratio' < 10 {
            local overfit_flag "Moderate concern"
        }
        if `heuristic_ratio' < 5 | `omitted' > 0 | `converged' != 1 {
            local overfit_flag "High concern"
        }

        post `post_age_diag' ("Decannulated-only home-any model") ("logit") ///
            (`decann_n') (`decann_home_n') (`decann_fac_n') (.) (.) (.) (.) ///
            (`k_slopes') (`omitted') (`converged') (`heuristic_ratio') (.) ///
            ("`overfit_flag'") ("ok")
        post `post_decann_notes' ("Interpret this subgroup descriptively only because the available information is limited.")
    }
}

post `post_age_notes' ("These age and discharge models are supplemental sensitivity analyses and should not be interpreted causally.")
post `post_age_notes' ("The home-any interaction model uses age_decade and averages predicted probabilities over the observed analytic cohort case mix.")
post `post_age_notes' ("`small_sample_note'")
post `post_age_notes' ("For binary models, support is summarized as the smaller of the event and nonevent counts divided by the number of free slope parameters.")
postclose `post_age_model'
postclose `post_age_diag'
postclose `post_age_notes'
postclose `post_decann_model'
postclose `post_decann_notes'

preserve
use `age_model_tbl', clear
keep if model_name == "Milestone-adjusted discharge ordinal model"
export excel using "`results_dir'/Supplemental Table - Age and Discharge Sensitivity Models.xlsx", ///
    sheet("Full-cohort ordered logit") firstrow(variables) replace
restore

preserve
use `age_model_tbl', clear
keep if model_name == "Age x decannulation home-any model"
export excel using "`results_dir'/Supplemental Table - Age and Discharge Sensitivity Models.xlsx", ///
    sheet("Interaction model") firstrow(variables) sheetreplace
restore

preserve
use `age_diag_tbl', clear
export excel using "`results_dir'/Supplemental Table - Age and Discharge Sensitivity Models.xlsx", ///
    sheet("Diagnostics") firstrow(variables) sheetreplace
restore

preserve
use `decann_desc_tbl', clear
export excel using "`results_dir'/Supplemental Table - Age and Discharge Sensitivity Models.xlsx", ///
    sheet("Decannulated subgroup") firstrow(variables) sheetreplace
restore

preserve
use `age_notes_tbl', clear
export excel using "`results_dir'/Supplemental Table - Age and Discharge Sensitivity Models.xlsx", ///
    sheet("Notes") firstrow(variables) sheetreplace
restore

preserve
use `decann_desc_tbl', clear
export excel using "`results_dir'/Supplemental Table - Decannulated Subgroup Descriptives.xlsx", ///
    sheet("Descriptives") firstrow(variables) replace
restore

preserve
use `decann_model_tbl', clear
export excel using "`results_dir'/Supplemental Table - Decannulated Subgroup Descriptives.xlsx", ///
    sheet("Model") firstrow(variables) sheetreplace
restore

preserve
use `decann_notes_tbl', clear
export excel using "`results_dir'/Supplemental Table - Decannulated Subgroup Descriptives.xlsx", ///
    sheet("Notes") firstrow(variables) sheetreplace
restore

preserve
tempvar discharge_jit
generate double `discharge_jit' = discharge_to + cond(decannulate == 1, 0.08, -0.08) + (runiform() - 0.5) * 0.16

twoway ///
    (scatter `discharge_jit' age if decannulate == 0, msymbol(Oh) mcolor(gs8) msize(medlarge)) ///
    (scatter `discharge_jit' age if decannulate == 1, msymbol(D) mcolor(black) msize(medlarge)), ///
    ylabel(1 "LTAC" 2 "SNF" 3 "Home w/ HH" 4 "Home", angle(0) labsize(medsmall)) ///
    yscale(range(0.6 4.4)) ///
    xtitle("Age (years)", size(medlarge)) ///
    ytitle("Discharge destination", size(medlarge)) ///
    title("Observed discharge disposition by age", size(large)) ///
    legend(order(1 "Not decannulated" 2 "Decannulated") rows(1) size(small) position(6)) ///
    note("`figure_note_age_raw'", size(vsmall)) ///
    scheme(s1mono) ///
    name(gr_age_raw, replace)
restore

if `interaction_fit_ok' {
    estimates restore age_home_interaction_model
    margins decannulate, at(age_decade = (2(0.5)8))
    marginsplot, ///
        recast(line) recastci(rarea) ///
        plot1opts(lcolor(black) lwidth(medthick) lpattern(solid)) ///
        plot2opts(lcolor(gs8) lwidth(medthick) lpattern(dash)) ///
        ci1opts(color(gs12)) ci2opts(color(gs14)) ///
        xlabel(2 "20" 2.5 "25" 3 "30" 3.5 "35" 4 "40" 4.5 "45" 5 "50" 5.5 "55" 6 "60" 6.5 "65" 7 "70" 7.5 "75" 8 "80", ///
            angle(45) labsize(vsmall)) ///
        xtitle("Age (years)", size(medlarge)) ///
        ytitle("Predicted Pr(home_any = 1)", size(medlarge)) ///
        title("Estimated probability of discharge to home", size(large)) ///
        subtitle("By decannulation status; averaged over observed case mix", size(small)) ///
        legend(order(1 "Not decannulated" 2 "Decannulated") rows(1) size(small) position(6)) ///
        note("`figure_note_age_model'", size(vsmall)) ///
        scheme(s1mono) ///
        name(gr_age_margins, replace)
}
else {
    twoway scatteri 0 0, ///
        msymbol(i) ///
        xlabel(none) ylabel(none) ///
        xscale(off) yscale(off) ///
        text(0.05 0 "Interaction model omitted", placement(c) size(medlarge)) ///
        text(-0.05 0 "See diagnostics worksheet and log.", placement(c) size(medlarge)) ///
        title("Estimated probability of discharge to home", size(large)) ///
        note("`figure_note_age_model'", size(vsmall)) ///
        legend(off) ///
        scheme(s1mono) ///
        name(gr_age_margins, replace)
}

graph combine gr_age_raw gr_age_margins, ///
    cols(2) xsize(12) ysize(5.5) ///
    imargin(small) ///
    name(gr_age_discharge, replace)

graph save gr_age_discharge ///
    "`results_dir'/Supplemental Figure - Age, Decannulation, and Discharge.gph", replace
quietly _export_graph_tiff, graphname(gr_age_discharge) ///
    outfile("`results_dir'/Supplemental Figure - Age, Decannulation, and Discharge.tiff") ///
    width(3600) height(1800)

di as txt "Supplemental sensitivity analyses completed."
di as txt "Expected outputs written to `results_dir'."

log close

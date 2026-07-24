* Canonical NRH-SCI-Vent orchestration entry point.

version 18
clear
set more off

args profile data_dir output_root run_id

if strtrim("`profile'") == "" local profile "smoke"
local profile = lower(strtrim("`profile'"))

if !inlist("`profile'", "smoke", "full", "release") {
    di as err "Invalid profile `profile'. Choose smoke, full, or release."
    exit 198
}

capture confirm file "PROJECT.yml"
if _rc {
    di as err "PROJECT.yml was not found. Run run_all.do from the repository root."
    exit 601
}

foreach required_file in ///
    "NRH SCI Cohort Preprocessing.do" ///
    "NRH SCI Cohort Paper Analysis.do" ///
    "NRH SCI Cohort Supplemental Sensitivity Analyses.do" ///
    "config/default.do" ///
    "config/smoke.do" ///
    "config/full.do" ///
    "config/release.do" ///
    "validation/expected_results_schema.csv" {
    capture confirm file "`required_file'"
    if _rc {
        di as err "Repository preflight could not find `required_file'."
        exit 601
    }
}

include "config/default.do"
if "`profile'" == "smoke" include "config/smoke.do"
if "`profile'" == "full" include "config/full.do"
if "`profile'" == "release" include "config/release.do"

if "`profile'" == "smoke" {
    if strtrim("`data_dir'") != "" & strtrim("`data_dir'") != "`nrh_smoke_data_dir'" {
        di as err "The smoke profile is restricted to the public synthetic data directory."
        exit 198
    }
    local data_dir "`nrh_smoke_data_dir'"
}
else {
    if strtrim("`data_dir'") == "" {
        di as err "The `profile' profile requires an explicit approved data directory."
        di as err "Restricted-data profiles never fall back to the synthetic fixture."
        exit 198
    }

    local nrh_data_dir_normalized = lower(subinstr(strtrim("`data_dir'"), "\", "/", .))
    while substr("`nrh_data_dir_normalized'", 1, 2) == "./" {
        local nrh_data_dir_normalized = substr("`nrh_data_dir_normalized'", 3, .)
    }
    while strpos("`nrh_data_dir_normalized'", "/./") {
        local nrh_data_dir_normalized = subinstr("`nrh_data_dir_normalized'", "/./", "/", .)
    }
    while regexm("`nrh_data_dir_normalized'", "/[.]$") {
        local nrh_data_dir_normalized = substr("`nrh_data_dir_normalized'", 1, ///
            strlen("`nrh_data_dir_normalized'") - 2)
    }
    local nrh_data_dir_normalized = regexr("`nrh_data_dir_normalized'", "/+$", "")
    if regexm("`nrh_data_dir_normalized'", "(^|/)[.][.]($|/)") {
        di as err "The `profile' profile data directory may not contain parent-directory segments."
        di as err "Pass an explicit normalized approved restricted-data directory."
        exit 198
    }
    if inlist("`nrh_data_dir_normalized'", "data/synthetic", "./data/synthetic") | ///
        regexm("`nrh_data_dir_normalized'", "/data/synthetic$") {
        di as err "The `profile' profile cannot use the public synthetic data directory."
        di as err "Pass an explicit approved restricted-data directory."
        exit 198
    }
}

if strtrim("`output_root'") == "" local output_root "`nrh_default_output_root'"

if strtrim("`run_id'") == "" {
    local nrh_date_num = date("`c(current_date)'", "DMY")
    local nrh_date : display %tdCCYY-NN-DD `nrh_date_num'
    local nrh_date = strtrim("`nrh_date'")
    local nrh_time = subinstr("`c(current_time)'", ":", "", .)
    local run_id "`nrh_date'T`nrh_time'_`profile'"
}

if strlen("`run_id'") > 128 | ///
    inlist("`run_id'", ".", "..") | ///
    strpos("`run_id'", "..") | ///
    !regexm("`run_id'", "^[A-Za-z0-9][A-Za-z0-9._-]*$") {
    di as err "run_id must be 1-128 path-safe characters: letters, numbers, dot, underscore, or hyphen."
    exit 198
}

local results_dir "`output_root'/`run_id'"
local log_dir "`results_dir'/Logs"
local master_log "`results_dir'/run_all.log"
local manifest_file "`results_dir'/run_manifest.csv"
local nrh_started "`c(current_date)' `c(current_time)'"

capture mkdir "`output_root'"
capture mkdir "`results_dir'"
if _rc {
    di as err "Could not create the new run directory `results_dir'."
    di as err "Run directories are never overwritten; choose a new run_id and a writable output root."
    exit 693
}

capture mkdir "`log_dir'"
if _rc {
    di as err "Could not create the run log directory."
    exit 693
}

tempname nrh_probe
capture file open `nrh_probe' using "`results_dir'/.nrh-write-test", write text replace
if _rc {
    di as err "The run directory is not writable."
    exit 603
}
file write `nrh_probe' "NRH-SCI-Vent write test" _n
file close `nrh_probe'
capture erase "`results_dir'/.nrh-write-test"

local nrh_overall_rc 0
local nrh_failed_stage ""
local nrh_preflight_rc .
local nrh_preprocessing_rc .
local nrh_cleaned_validation_rc .
local nrh_paper_rc .
local nrh_supplemental_rc .
local nrh_output_check_rc .
local nrh_manifest_rc .
local nrh_toplog_rc 0

* Write the top-level log explicitly so it cannot echo private paths or results.
tempname nrh_toplog
capture file open `nrh_toplog' using "`master_log'", write text replace
local nrh_toplog_rc = _rc
if `nrh_toplog_rc' {
    di as err "Could not create the value-free top-level run log."
    exit `nrh_toplog_rc'
}
file write `nrh_toplog' "NRH-SCI-Vent orchestration" _n
file write `nrh_toplog' "profile=`profile'" _n
file write `nrh_toplog' "run_id=`run_id'" _n
file write `nrh_toplog' "started_local=`nrh_started'" _n
file close `nrh_toplog'

di as txt "NRH-SCI-Vent orchestration"
di as txt "Profile: `profile'"
di as txt "Run ID: `run_id'"
di as txt "Started: `nrh_started'"

* Stage 1: profile and input preflight.
if `nrh_overall_rc' == 0 {
    local nrh_preflight_rc 0

    capture confirm file "`data_dir'/`nrh_source_file'"
    if _rc {
        local nrh_preflight_rc 601
        if "`profile'" == "smoke" {
            di as err "The public synthetic fixture is not present at `data_dir'/`nrh_source_file'."
            di as err "NRH-005 will supply the approved synthetic fixture; no restricted-data fallback was attempted."
        }
        else {
            di as err "Could not find the required source CSV in the approved data directory."
        }
    }

    if `nrh_preflight_rc' == 0 & `nrh_refuse_derived_overwrite' {
        capture confirm new file "`data_dir'/nrh-sci-raw.dta"
        if _rc local nrh_preflight_rc 602
        capture confirm new file "`data_dir'/nrh-sci-cleaned.dta"
        if _rc local nrh_preflight_rc 602
        if `nrh_preflight_rc' {
            di as err "Release mode refuses to overwrite existing derived Stata datasets."
            di as err "Use a new approved data-staging directory."
        }
    }

    di as txt "NRH_STAGE preflight rc=`nrh_preflight_rc'"
    capture file open `nrh_toplog' using "`master_log'", write text append
    if _rc {
        local nrh_toplog_rc = _rc
    }
    else {
        file write `nrh_toplog' "stage.preflight.rc=`nrh_preflight_rc'" _n
        file close `nrh_toplog'
    }
    if `nrh_preflight_rc' {
        local nrh_overall_rc `nrh_preflight_rc'
        local nrh_failed_stage "preflight"
    }
}

* Stage 2: unchanged preprocessing body through the legacy entry point.
if `nrh_overall_rc' == 0 {
    capture noisily do "NRH SCI Cohort Preprocessing.do" ///
        "`data_dir'" "`output_root'" "`run_id'"
    local nrh_preprocessing_rc = _rc
    capture log close
    di as txt "NRH_STAGE preprocessing rc=`nrh_preprocessing_rc'"
    capture file open `nrh_toplog' using "`master_log'", write text append
    if _rc {
        local nrh_toplog_rc = _rc
    }
    else {
        file write `nrh_toplog' "stage.preprocessing.rc=`nrh_preprocessing_rc'" _n
        file close `nrh_toplog'
    }
    if `nrh_preprocessing_rc' {
        local nrh_overall_rc `nrh_preprocessing_rc'
        local nrh_failed_stage "preprocessing"
    }
}

* Stage 3: sentinel cleaned-file validation. NRH-003/006 add full contracts.
if `nrh_overall_rc' == 0 {
    capture confirm file "`data_dir'/nrh-sci-cleaned.dta"
    local nrh_cleaned_validation_rc = _rc
    di as txt "NRH_STAGE cleaned_validation rc=`nrh_cleaned_validation_rc'"
    capture file open `nrh_toplog' using "`master_log'", write text append
    if _rc {
        local nrh_toplog_rc = _rc
    }
    else {
        file write `nrh_toplog' "stage.cleaned_validation.rc=`nrh_cleaned_validation_rc'" _n
        file close `nrh_toplog'
    }
    if `nrh_cleaned_validation_rc' {
        local nrh_overall_rc `nrh_cleaned_validation_rc'
        local nrh_failed_stage "cleaned_validation"
    }
}

* Stage 4: unchanged paper-analysis body through the legacy entry point.
if `nrh_overall_rc' == 0 {
    capture noisily do "NRH SCI Cohort Paper Analysis.do" ///
        "`data_dir'" "`output_root'" "`run_id'"
    local nrh_paper_rc = _rc
    capture log close
    di as txt "NRH_STAGE paper rc=`nrh_paper_rc'"
    capture file open `nrh_toplog' using "`master_log'", write text append
    if _rc {
        local nrh_toplog_rc = _rc
    }
    else {
        file write `nrh_toplog' "stage.paper.rc=`nrh_paper_rc'" _n
        file close `nrh_toplog'
    }
    if `nrh_paper_rc' {
        local nrh_overall_rc `nrh_paper_rc'
        local nrh_failed_stage "paper"
    }
}

* Stage 5: unchanged supplemental-analysis body through the legacy entry point.
if `nrh_overall_rc' == 0 {
    local nrh_plus_before "`c(sysdir_plus)'"
    capture noisily do "NRH SCI Cohort Supplemental Sensitivity Analyses.do" ///
        "`data_dir'" "`output_root'" "`run_id'"
    local nrh_supplemental_rc = _rc
    capture log close
    quietly sysdir set PLUS "`nrh_plus_before'"
    capture adopath - "`c(tmpdir)'codex_oparallel_plus"
    di as txt "NRH_STAGE supplemental rc=`nrh_supplemental_rc'"
    capture file open `nrh_toplog' using "`master_log'", write text append
    if _rc {
        local nrh_toplog_rc = _rc
    }
    else {
        file write `nrh_toplog' "stage.supplemental.rc=`nrh_supplemental_rc'" _n
        file close `nrh_toplog'
    }
    if `nrh_supplemental_rc' {
        local nrh_overall_rc `nrh_supplemental_rc'
        local nrh_failed_stage "supplemental"
    }
}

* Stage 6: check the output-presence rows in the public value-free contract.
if `nrh_overall_rc' == 0 {
    local nrh_output_check_rc 0
    capture noisily import delimited using "`nrh_expected_results_contract'", ///
        clear varnames(1) stringcols(_all) bindquote(strict)
    if _rc {
        local nrh_output_check_rc = _rc
    }
    else {
        capture confirm variable result_family term_or_level required
        if _rc {
            local nrh_output_check_rc = _rc
            di as err "The expected-results contract is missing required output-inventory columns."
        }
        else {
            quietly keep if result_family == "output_presence" & required == "true"
            quietly count
            if r(N) == 0 {
                local nrh_output_check_rc 459
            }
            else {
                local nrh_output_count = r(N)
                forvalues nrh_i = 1/`nrh_output_count' {
                    local nrh_expected_file = term_or_level[`nrh_i']
                    capture confirm file "`results_dir'/`nrh_expected_file'"
                    if _rc {
                        di as err "Expected output is missing: `nrh_expected_file'"
                        local nrh_output_check_rc 601
                    }
                }
            }
        }
    }

    di as txt "NRH_STAGE output_check rc=`nrh_output_check_rc'"
    capture file open `nrh_toplog' using "`master_log'", write text append
    if _rc {
        local nrh_toplog_rc = _rc
    }
    else {
        file write `nrh_toplog' "stage.output_check.rc=`nrh_output_check_rc'" _n
        file close `nrh_toplog'
    }
    if `nrh_output_check_rc' {
        local nrh_overall_rc `nrh_output_check_rc'
        local nrh_failed_stage "output_check"
    }
}

if `nrh_toplog_rc' & `nrh_overall_rc' == 0 {
    local nrh_overall_rc `nrh_toplog_rc'
    local nrh_failed_stage "top_log"
}

* Stage 7: always attempt a value-free runtime manifest after run creation.
local nrh_completed "`c(current_date)' `c(current_time)'"
tempname nrh_manifest
capture file open `nrh_manifest' using "`manifest_file'", write text replace
local nrh_manifest_rc = _rc
if `nrh_manifest_rc' == 0 {
    file write `nrh_manifest' "key,value" _n
    file write `nrh_manifest' "manifest_version,1" _n
    file write `nrh_manifest' "project,NRH-SCI-Vent" _n
    file write `nrh_manifest' "profile,`profile'" _n
    file write `nrh_manifest' "run_id,`run_id'" _n
    file write `nrh_manifest' "started_local,`nrh_started'" _n
    file write `nrh_manifest' "completed_local,`nrh_completed'" _n
    file write `nrh_manifest' "stata_version,`c(stata_version)'" _n
    file write `nrh_manifest' "stata_edition,`c(edition_real)'" _n
    file write `nrh_manifest' "stata_build_date,`c(born_date)'" _n
    file write `nrh_manifest' "operating_system,`c(os)'" _n
    file write `nrh_manifest' "preflight_rc,`nrh_preflight_rc'" _n
    file write `nrh_manifest' "preprocessing_rc,`nrh_preprocessing_rc'" _n
    file write `nrh_manifest' "cleaned_validation_rc,`nrh_cleaned_validation_rc'" _n
    file write `nrh_manifest' "paper_rc,`nrh_paper_rc'" _n
    file write `nrh_manifest' "supplemental_rc,`nrh_supplemental_rc'" _n
    file write `nrh_manifest' "output_check_rc,`nrh_output_check_rc'" _n
    file write `nrh_manifest' "top_log_rc,`nrh_toplog_rc'" _n
    file write `nrh_manifest' "failed_stage,`nrh_failed_stage'" _n
    file write `nrh_manifest' "overall_rc,`nrh_overall_rc'" _n
    file close `nrh_manifest'
}
else if `nrh_overall_rc' == 0 {
    local nrh_overall_rc `nrh_manifest_rc'
    local nrh_failed_stage "manifest"
}

di as txt "NRH_STAGE manifest rc=`nrh_manifest_rc'"
di as txt "NRH_RUN_COMPLETE rc=`nrh_overall_rc' failed_stage=`nrh_failed_stage'"
capture file open `nrh_toplog' using "`master_log'", write text append
if !_rc {
    file write `nrh_toplog' "stage.manifest.rc=`nrh_manifest_rc'" _n
    file write `nrh_toplog' "completed_local=`nrh_completed'" _n
    file write `nrh_toplog' "overall.rc=`nrh_overall_rc'" _n
    file write `nrh_toplog' "failed_stage=`nrh_failed_stage'" _n
    file close `nrh_toplog'
}

exit `nrh_overall_rc'

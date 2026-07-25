* Canonical NRH-SCI-Vent orchestration entry point.

version 18
clear
set more off
global NRH_DEPENDENCY_MANIFEST_SHA256 ""

args profile data_dir output_root run_id

if strtrim("`profile'") == "" local profile "verify"
local profile = lower(strtrim("`profile'"))

if "`profile'" == "smoke" {
    di as err "The public smoke profile is retired because no synthetic cohort is authorized."
    di as err "Use the verify profile for explicitly no-data public verification."
    exit 198
}

if !inlist("`profile'", "verify", "full", "release") {
    di as err "Invalid profile `profile'. Choose verify, full, or release."
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
    "config/verify.do" ///
    "config/full.do" ///
    "config/release.do" ///
    "code/00_preflight.do" ///
    "code/lib/data_contracts.do" ///
    "code/lib/value_mappings.do" ///
    "tests/test_public_contracts.do" ///
    "tests/test_value_mappings.do" ///
    "vendor/stata/manifest.csv" ///
    "validation/data_sources.csv" ///
    "validation/expected_results_schema.csv" ///
    "validation/sample_flow_schema.csv" ///
    "validation/source_schema.csv" ///
    "validation/value_mappings.csv" ///
    "validation/validation_rules.csv" {
    capture confirm file "`required_file'"
    if _rc {
        di as err "Repository preflight could not find `required_file'."
        exit 601
    }
}

include "config/default.do"
if "`profile'" == "verify" include "config/verify.do"
if "`profile'" == "full" include "config/full.do"
if "`profile'" == "release" include "config/release.do"

local nrh_run_scope "restricted_analysis"
local nrh_data_accessed "false"

if "`profile'" == "verify" {
    local nrh_run_scope "no_data"
    if strtrim("`data_dir'") != "" {
        di as err "The verify profile does not accept or access a data directory."
        exit 198
    }
}
else {
    if strtrim("`data_dir'") == "" {
        di as err "The `profile' profile requires an explicit approved data directory."
        di as err "Restricted-data profiles never infer or substitute another input."
        exit 198
    }

    local nrh_data_dir_normalized = lower(subinstr(strtrim("`data_dir'"), "\", "/", .))
    while strpos("`nrh_data_dir_normalized'", "//") {
        local nrh_data_dir_normalized = subinstr("`nrh_data_dir_normalized'", "//", "/", .)
    }
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
        di as err "The retired public-fixture location is not an authorized clinical-data directory."
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
local nrh_dependency_preflight_rc .
local nrh_dependency_manifest_sha256 ""
local nrh_public_contract_rc .
local nrh_public_contract_status "not_run"
local nrh_mapping_unit_rc .
local nrh_mapping_unit_status "not_run"
local nrh_preflight_rc .
local nrh_preflight_status "not_run"
local nrh_source_contract_rc .
local nrh_source_contract_version .
local nrh_source_contract_status "not_run"
local nrh_preprocessing_rc .
local nrh_preprocessing_status "not_run"
local nrh_cleaned_validation_rc .
local nrh_cleaned_validation_status "not_run"
local nrh_paper_rc .
local nrh_paper_status "not_run"
local nrh_supplemental_rc .
local nrh_supplemental_status "not_run"
local nrh_output_check_rc .
local nrh_output_check_status "not_run"
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
file write `nrh_toplog' "run_scope=`nrh_run_scope'" _n
file write `nrh_toplog' "run_id=`run_id'" _n
file write `nrh_toplog' "started_local=`nrh_started'" _n
file close `nrh_toplog'

di as txt "NRH-SCI-Vent orchestration"
di as txt "Profile: `profile'"
di as txt "Run ID: `run_id'"
di as txt "Started: `nrh_started'"

* Stage 1: verify the controlled Stata environment before any data access.
if `nrh_overall_rc' == 0 {
    capture noisily do "code/00_preflight.do" ///
        "`log_dir'/dependency_preflight.log"
    local nrh_dependency_preflight_rc = _rc
    local nrh_dependency_manifest_sha256 "$NRH_DEPENDENCY_MANIFEST_SHA256"
    di as txt "NRH_STAGE dependency_preflight rc=`nrh_dependency_preflight_rc'"
    capture file open `nrh_toplog' using "`master_log'", write text append
    if _rc {
        local nrh_toplog_rc = _rc
    }
    else {
        file write `nrh_toplog' ///
            "stage.dependency_preflight.rc=`nrh_dependency_preflight_rc'" _n
        file close `nrh_toplog'
    }
    if `nrh_dependency_preflight_rc' {
        local nrh_overall_rc `nrh_dependency_preflight_rc'
        local nrh_failed_stage "dependency_preflight"
    }
}

* Stages 2-3 for verify: public contracts and atomic mapping behavior only.
if "`profile'" == "verify" & `nrh_overall_rc' == 0 {
    capture noisily do "tests/test_public_contracts.do"
    local nrh_public_contract_rc = _rc
    local nrh_public_contract_status = ///
        cond(`nrh_public_contract_rc' == 0, "passed", "failed")
    di as txt "NRH_STAGE public_contract rc=`nrh_public_contract_rc'"
    capture file open `nrh_toplog' using "`master_log'", write text append
    if _rc {
        local nrh_toplog_rc = _rc
    }
    else {
        file write `nrh_toplog' ///
            "stage.public_contract.rc=`nrh_public_contract_rc'" _n
        file close `nrh_toplog'
    }
    if `nrh_public_contract_rc' {
        local nrh_overall_rc `nrh_public_contract_rc'
        local nrh_failed_stage "public_contract"
    }
}

if "`profile'" == "verify" & `nrh_overall_rc' == 0 {
    capture noisily do "tests/test_value_mappings.do"
    local nrh_mapping_unit_rc = _rc
    local nrh_mapping_unit_status = ///
        cond(`nrh_mapping_unit_rc' == 0, "passed", "failed")
    di as txt "NRH_STAGE mapping_unit rc=`nrh_mapping_unit_rc'"
    capture file open `nrh_toplog' using "`master_log'", write text append
    if _rc {
        local nrh_toplog_rc = _rc
    }
    else {
        file write `nrh_toplog' ///
            "stage.mapping_unit.rc=`nrh_mapping_unit_rc'" _n
        file close `nrh_toplog'
    }
    if `nrh_mapping_unit_rc' {
        local nrh_overall_rc `nrh_mapping_unit_rc'
        local nrh_failed_stage "mapping_unit"
    }
}

* Data stage 1: profile and input preflight.
if "`profile'" != "verify" & `nrh_overall_rc' == 0 {
    local nrh_preflight_rc 0

    capture confirm file "`data_dir'/`nrh_source_file'"
    if _rc {
        local nrh_preflight_rc 601
        di as err "Could not find the required source CSV in the approved data directory."
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

    local nrh_preflight_status = ///
        cond(`nrh_preflight_rc' == 0, "passed", "failed")
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

* Data stage 2: validate the ordered source contract before preprocessing.
if "`profile'" != "verify" & `nrh_overall_rc' == 0 {
    capture noisily do "code/lib/data_contracts.do"
    local nrh_source_contract_rc = _rc
    if `nrh_source_contract_rc' == 0 {
        local nrh_data_accessed "true"
        capture noisily nrh_validate_source_contract ///
            using "`data_dir'/`nrh_source_file'", ///
            schema("validation/source_schema.csv") ///
            rules("validation/validation_rules.csv") ///
            mode("strict") ///
            log("`log_dir'/source_contract.log")
        local nrh_source_contract_rc = _rc
        if `nrh_source_contract_rc' == 0 {
            local nrh_source_contract_version = r(contract_version)
        }
    }
    local nrh_source_contract_status = ///
        cond(`nrh_source_contract_rc' == 0, "passed", "failed")
    di as txt "NRH_STAGE source_contract rc=`nrh_source_contract_rc'"
    capture file open `nrh_toplog' using "`master_log'", write text append
    if _rc {
        local nrh_toplog_rc = _rc
    }
    else {
        file write `nrh_toplog' ///
            "stage.source_contract.rc=`nrh_source_contract_rc'" _n
        file close `nrh_toplog'
    }
    if `nrh_source_contract_rc' {
        local nrh_overall_rc `nrh_source_contract_rc'
        local nrh_failed_stage "source_contract"
    }
}

* Data stage 3: preprocessing through the legacy entry point.
if "`profile'" != "verify" & `nrh_overall_rc' == 0 {
    capture noisily do "NRH SCI Cohort Preprocessing.do" ///
        "`data_dir'" "`output_root'" "`run_id'"
    local nrh_preprocessing_rc = _rc
    local nrh_preprocessing_status = ///
        cond(`nrh_preprocessing_rc' == 0, "passed", "failed")
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

* Data stage 4: sentinel cleaned-file validation. NRH-006 adds the contract.
if "`profile'" != "verify" & `nrh_overall_rc' == 0 {
    capture confirm file "`data_dir'/nrh-sci-cleaned.dta"
    local nrh_cleaned_validation_rc = _rc
    local nrh_cleaned_validation_status = ///
        cond(`nrh_cleaned_validation_rc' == 0, "passed", "failed")
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

* Data stage 5: paper analysis through the legacy entry point.
if "`profile'" != "verify" & `nrh_overall_rc' == 0 {
    capture noisily do "NRH SCI Cohort Paper Analysis.do" ///
        "`data_dir'" "`output_root'" "`run_id'"
    local nrh_paper_rc = _rc
    local nrh_paper_status = cond(`nrh_paper_rc' == 0, "passed", "failed")
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

* Data stage 6: supplemental analysis through the legacy entry point.
if "`profile'" != "verify" & `nrh_overall_rc' == 0 {
    capture noisily do "NRH SCI Cohort Supplemental Sensitivity Analyses.do" ///
        "`data_dir'" "`output_root'" "`run_id'"
    local nrh_supplemental_rc = _rc
    local nrh_supplemental_status = ///
        cond(`nrh_supplemental_rc' == 0, "passed", "failed")
    capture log close
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

* Data stage 7: check required output presence.
if "`profile'" != "verify" & `nrh_overall_rc' == 0 {
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

    local nrh_output_check_status = ///
        cond(`nrh_output_check_rc' == 0, "passed", "failed")
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

* Stage 9: always attempt a value-free runtime manifest after run creation.
local nrh_completed "`c(current_date)' `c(current_time)'"
tempname nrh_manifest
capture file open `nrh_manifest' using "`manifest_file'", write text replace
local nrh_manifest_rc = _rc
if `nrh_manifest_rc' == 0 {
    file write `nrh_manifest' "key,value" _n
    file write `nrh_manifest' "manifest_version,2" _n
    file write `nrh_manifest' "project,NRH-SCI-Vent" _n
    file write `nrh_manifest' "profile,`profile'" _n
    file write `nrh_manifest' "run_scope,`nrh_run_scope'" _n
    file write `nrh_manifest' "data_accessed,`nrh_data_accessed'" _n
    file write `nrh_manifest' "run_id,`run_id'" _n
    file write `nrh_manifest' "started_local,`nrh_started'" _n
    file write `nrh_manifest' "completed_local,`nrh_completed'" _n
    file write `nrh_manifest' "stata_version,`c(stata_version)'" _n
    file write `nrh_manifest' "stata_edition,`c(edition_real)'" _n
    file write `nrh_manifest' "stata_build_date,`c(born_date)'" _n
    file write `nrh_manifest' "operating_system,`c(os)'" _n
    file write `nrh_manifest' ///
        "dependency_manifest_sha256,`nrh_dependency_manifest_sha256'" _n
    file write `nrh_manifest' "dependency_vendor_path,vendor/stata/plus" _n
    file write `nrh_manifest' ///
        "dependency_preflight_rc,`nrh_dependency_preflight_rc'" _n
    file write `nrh_manifest' ///
        "public_contract_rc,`nrh_public_contract_rc'" _n
    file write `nrh_manifest' ///
        "public_contract_status,`nrh_public_contract_status'" _n
    file write `nrh_manifest' "mapping_unit_rc,`nrh_mapping_unit_rc'" _n
    file write `nrh_manifest' ///
        "mapping_unit_status,`nrh_mapping_unit_status'" _n
    file write `nrh_manifest' "preflight_rc,`nrh_preflight_rc'" _n
    file write `nrh_manifest' "preflight_status,`nrh_preflight_status'" _n
    file write `nrh_manifest' ///
        "source_contract_version,`nrh_source_contract_version'" _n
    file write `nrh_manifest' ///
        "source_contract_rc,`nrh_source_contract_rc'" _n
    file write `nrh_manifest' ///
        "source_contract_status,`nrh_source_contract_status'" _n
    file write `nrh_manifest' "preprocessing_rc,`nrh_preprocessing_rc'" _n
    file write `nrh_manifest' ///
        "preprocessing_status,`nrh_preprocessing_status'" _n
    file write `nrh_manifest' "cleaned_validation_rc,`nrh_cleaned_validation_rc'" _n
    file write `nrh_manifest' ///
        "cleaned_validation_status,`nrh_cleaned_validation_status'" _n
    file write `nrh_manifest' "paper_rc,`nrh_paper_rc'" _n
    file write `nrh_manifest' "paper_status,`nrh_paper_status'" _n
    file write `nrh_manifest' "supplemental_rc,`nrh_supplemental_rc'" _n
    file write `nrh_manifest' ///
        "supplemental_status,`nrh_supplemental_status'" _n
    file write `nrh_manifest' "output_check_rc,`nrh_output_check_rc'" _n
    file write `nrh_manifest' ///
        "output_check_status,`nrh_output_check_status'" _n
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
    file write `nrh_toplog' "data_accessed=`nrh_data_accessed'" _n
    file write `nrh_toplog' ///
        "stage.public_contract.status=`nrh_public_contract_status'" _n
    file write `nrh_toplog' ///
        "stage.mapping_unit.status=`nrh_mapping_unit_status'" _n
    file write `nrh_toplog' ///
        "stage.preflight.status=`nrh_preflight_status'" _n
    file write `nrh_toplog' ///
        "stage.source_contract.status=`nrh_source_contract_status'" _n
    file write `nrh_toplog' ///
        "stage.preprocessing.status=`nrh_preprocessing_status'" _n
    file write `nrh_toplog' ///
        "stage.cleaned_validation.status=`nrh_cleaned_validation_status'" _n
    file write `nrh_toplog' "stage.paper.status=`nrh_paper_status'" _n
    file write `nrh_toplog' ///
        "stage.supplemental.status=`nrh_supplemental_status'" _n
    file write `nrh_toplog' ///
        "stage.output_check.status=`nrh_output_check_status'" _n
    file write `nrh_toplog' "overall.rc=`nrh_overall_rc'" _n
    file write `nrh_toplog' "failed_stage=`nrh_failed_stage'" _n
    file close `nrh_toplog'
}

capture macro drop NRH_DEPENDENCY_MANIFEST_SHA256
exit `nrh_overall_rc'

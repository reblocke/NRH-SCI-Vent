* Data-free structural tests for the public source and baseline contracts.

version 18
clear
set more off

args repo_root status_file
if strtrim("`repo_root'") == "" local repo_root "`c(pwd)'"
quietly cd "`repo_root'"

* Both executable contract libraries must load without opening source data.
do "code/lib/data_contracts.do"
capture program list nrh_validate_source_contract
assert _rc == 0

do "code/lib/value_mappings.do"
capture program list nrh_apply_value_mapping
assert _rc == 0

tempfile nrh_schema_fields

quietly import delimited using "validation/source_schema.csv", clear ///
    varnames(1) stringcols(_all) bindquote(strict) encoding("UTF-8")
confirm variable source_field required storage_type semantic_type ///
    allowed_values missing_allowed sensitivity used_by notes
assert _N == 36
isid source_field
assert required == "true"
assert inlist(storage_type, "numeric", "string")
assert inlist(missing_allowed, "true", "false")
assert inlist(semantic_type, "identifier", "integer", "date_mdy", ///
    "categorical", "numeric_with_exceptions", "free_text", ///
    "empty_placeholder")
assert inlist(sensitivity, "direct_identifier", "quasi_identifier", ///
    "sensitive_clinical", "restricted_free_text", "operational")
keep source_field
save `nrh_schema_fields', replace

quietly import delimited using "validation/validation_rules.csv", clear ///
    varnames(1) stringcols(_all) bindquote(strict) encoding("UTF-8")
confirm variable rule_id source_field rule_type parameters strict_action ///
    development_action disclosure_output notes
isid rule_id
assert inlist(rule_type, "unique_nonmissing", "parseable_date", "integer", ///
    "nonnegative", "all_missing", "numeric_or_allowed_literal")
assert strict_action == "fail"
assert inlist(development_action, "fail", "warn")
assert disclosure_output == "aggregate_count_only"
merge m:1 source_field using `nrh_schema_fields'
assert _merge != 1
drop _merge

quietly import delimited using "validation/data_sources.csv", clear ///
    varnames(1) stringcols(_all) bindquote(strict) encoding("UTF-8")
confirm variable source_id file_name format delimiter encoding ///
    line_endings header_row expected_columns classification ///
    contract_version approval_status approval_date approver_role ///
    secure_review_artifact_id notes
assert _N == 1
isid source_id
assert file_name == "Working NRH SCI.csv"
assert classification == "restricted clinical data"
assert contract_version == "2"
assert approval_status == "approved"

quietly import delimited using "validation/value_mappings.csv", clear ///
    varnames(1) stringcols(_all) bindquote(strict) encoding("UTF-8")
confirm variable mapping_id source_field input_kind normalized_literal ///
    target_variable action target_code target_label decision_id ///
    unlisted_action notes
isid mapping_id
assert inlist(input_kind, "literal", "blank")
assert inlist(action, "assign", "set_missing")
assert unlisted_action == "reject"
assert normalized_literal == "" if input_kind == "blank"
assert action == "set_missing" if input_kind == "blank"
assert target_code == "" if action == "set_missing"
assert target_label == "" if action == "set_missing"
merge m:1 source_field using `nrh_schema_fields'
assert _merge != 1
drop _merge

quietly import delimited using "validation/expected_results_schema.csv", ///
    clear varnames(1) stringcols(_all) bindquote(strict) encoding("UTF-8")
confirm variable result_id expected_value comparison_method ///
    absolute_tolerance relative_tolerance required disclosure_class
isid result_id
assert strtrim(expected_value) == ""

quietly import delimited using "validation/sample_flow_schema.csv", clear ///
    varnames(1) stringcols(_all) bindquote(strict) encoding("UTF-8")
confirm variable flow_id expected_before_n expected_excluded_n ///
    expected_after_n comparison_method absolute_tolerance required ///
    disclosure_class
isid flow_id
assert strtrim(expected_before_n) == ""
assert strtrim(expected_excluded_n) == ""
assert strtrim(expected_after_n) == ""

if strtrim("`status_file'") != "" {
    tempname nrh_status
    file open `nrh_status' using "`status_file'", write text replace
    file write `nrh_status' "0" _n
    file close `nrh_status'
}

di as result "Public contract structural tests passed without source data."

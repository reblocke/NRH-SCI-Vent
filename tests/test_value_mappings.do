* Synthetic-only tests for the NRH-004 value-mapping contract.

version 18
clear
set more off

args repo_root status_file
if strtrim("`repo_root'") == "" local repo_root "`c(pwd)'"
quietly cd "`repo_root'"

local nrh_contract "validation/value_mappings.csv"
capture confirm file "`nrh_contract'"
if _rc {
    di as err "The public value-mapping contract is unavailable."
    exit 601
}

do "code/lib/value_mappings.do"

tempfile nrh_contract_rows nrh_diagnostics
quietly import delimited using "`nrh_contract'", clear varnames(1) ///
    stringcols(_all) bindquote(strict) encoding("UTF-8")
confirm variable mapping_id source_field input_kind normalized_literal ///
    target_variable action target_code target_label decision_id ///
    unlisted_action notes

quietly egen long nrh_pair = group(source_field target_variable)
quietly summarize nrh_pair, meanonly
local nrh_pair_count = r(max)
save `nrh_contract_rows', replace

* Every approved literal and missing action must execute to its declared result.
forvalues nrh_pair = 1/`nrh_pair_count' {
    use `nrh_contract_rows', clear
    quietly keep if nrh_pair == `nrh_pair'
    local nrh_source_field `"`=source_field[1]'"'
    local nrh_target `"`=target_variable[1]'"'

    generate strL nrh_test_input = normalized_literal
    replace nrh_test_input = "" if input_kind == "blank"
    generate double nrh_expected_code = real(target_code)

    nrh_apply_value_mapping nrh_test_input, ///
        sourcefield("`nrh_source_field'") target(`nrh_target') ///
        contract("`nrh_contract'") diagnostics("`nrh_diagnostics'")

    assert `nrh_target' == nrh_expected_code if action == "assign"
    assert missing(`nrh_target') if action == "set_missing"

    decode `nrh_target', generate(nrh_actual_label)
    assert nrh_actual_label == target_label if action == "assign"
    assert nrh_actual_label == "" if action == "set_missing"
}

* The prospective compound domain is complete for all C1-C8 x AIS A-D cells
* and for both approved level-specific exception forms.
use `nrh_contract_rows', clear
forvalues nrh_level = 1/8 {
    foreach nrh_grade in A B C D {
        local nrh_grade_code = ///
            cond("`nrh_grade'" == "A", 1, ///
            cond("`nrh_grade'" == "B", 2, ///
            cond("`nrh_grade'" == "C", 3, 4)))

        quietly count if source_field == "asiaclassificationatrehab" & ///
            normalized_literal == "C`nrh_level' AIS `nrh_grade'" & ///
            target_variable == "level" & action == "assign" & ///
            real(target_code) == `nrh_level'
        assert r(N) == 1

        quietly count if source_field == "asiaclassificationatrehab" & ///
            normalized_literal == "C`nrh_level' AIS `nrh_grade'" & ///
            target_variable == "asia_class" & action == "assign" & ///
            real(target_code) == `nrh_grade_code'
        assert r(N) == 1
    }

    quietly count if source_field == "asiaclassificationatrehab" & ///
        normalized_literal == "C`nrh_level' AIA A" & ///
        target_variable == "level" & action == "assign" & ///
        real(target_code) == `nrh_level'
    assert r(N) == 1
    quietly count if source_field == "asiaclassificationatrehab" & ///
        normalized_literal == "C`nrh_level' AIA A" & ///
        target_variable == "asia_class" & action == "assign" & ///
        real(target_code) == 1
    assert r(N) == 1

    quietly count if source_field == "asiaclassificationatrehab" & ///
        normalized_literal == "C`nrh_level' CENTRAL CORD" & ///
        target_variable == "level" & action == "assign" & ///
        real(target_code) == `nrh_level'
    assert r(N) == 1
    quietly count if source_field == "asiaclassificationatrehab" & ///
        normalized_literal == "C`nrh_level' CENTRAL CORD" & ///
        target_variable == "asia_class" & action == "assign" & ///
        real(target_code) == 4
    assert r(N) == 1
}

quietly count if source_field == ///
    "weanoffventduringthedayandcontin" & ///
    target_variable == "partial_wean_at_admit" & ///
    input_kind == "blank" & action == "set_missing"
assert r(N) == 1
foreach nrh_unknown in UNKNOWN UNKOWN {
    quietly count if source_field == ///
        "weanoffventduringthedayandcontin" & ///
        target_variable == "partial_wean_at_admit" & ///
        normalized_literal == "`nrh_unknown'" & action == "set_missing"
    assert r(N) == 1
}

* Normalization accepts case and surrounding-space variants without changing
* the approved literal contract.
clear
input str20 nrh_test_input
" f "
"male"
"   "
end
nrh_apply_value_mapping nrh_test_input, sourcefield("sex") target(sex) ///
    contract("`nrh_contract'") diagnostics("`nrh_diagnostics'")
assert sex == 0 in 1
assert sex == 1 in 2
assert missing(sex) in 3
local nrh_sex_label_0 : label (sex) 0
local nrh_sex_label_1 : label (sex) 1
assert "`nrh_sex_label_0'" == "Female"
assert "`nrh_sex_label_1'" == "Male"

* A generated synthetic unlisted literal must fail closed without appearing in
* the aggregate diagnostics. The literal is assembled at runtime so it is not
* embedded in this public test source.
clear
set obs 1
generate strL nrh_test_input = ///
    char(85) + char(78) + char(76) + char(73) + char(83) + char(84) + ///
    char(69) + char(68)
local nrh_rejected_literal = nrh_test_input[1]
capture noisily nrh_apply_value_mapping nrh_test_input, ///
    sourcefield("sex") target(sex) contract("`nrh_contract'") ///
    diagnostics("`nrh_diagnostics'")
assert _rc == 459

tempname nrh_handle
local nrh_posix_users = "/" + "users" + "/"
local nrh_windows_users = char(92) + "users" + char(92)
file open `nrh_handle' using "`nrh_diagnostics'", read text
file read `nrh_handle' nrh_line
while r(eof) == 0 {
    assert strpos("`nrh_line'", "`nrh_rejected_literal'") == 0
    assert strpos(lower("`nrh_line'"), "row_number") == 0
    assert strpos(lower("`nrh_line'"), "`nrh_posix_users'") == 0
    assert strpos(lower("`nrh_line'"), "`nrh_windows_users'") == 0
    file read `nrh_handle' nrh_line
}
file close `nrh_handle'

if strtrim("`status_file'") != "" {
    tempname nrh_status
    file open `nrh_status' using "`status_file'", write text replace
    file write `nrh_status' "0" _n
    file close `nrh_status'
}

di as result "NRH value-mapping synthetic tests passed."

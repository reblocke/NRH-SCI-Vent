* Synthetic-only tests for the NRH-004 source contract version 2.

version 18
clear
set more off

args repo_root status_file
if strtrim("`repo_root'") == "" local repo_root "`c(pwd)'"
cd "`repo_root'"

capture program drop nrh_test_write_fixture
program define nrh_test_write_fixture
    version 18
    syntax using/, VARIANT(string)

    clear
    set obs 2

    generate long pat_mrn_id = _n
    generate byte age = 49 + _n
    generate str8 sex = "M"
    generate str16 dateofinjury = "01/01/2020"
    generate str8 injurylevelreportedpriortorehab = "C4"
    replace injurylevelreportedpriortorehab = "C8" in 2
    generate str8 ribfracturesyorn = "N"
    generate str8 pneumothoraxyorn = "N"
    generate str40 chesttubeyorn = "N"
    generate byte daysfrominjurytointubation = 0
    generate str20 reasonforintubation = "TEST"
    generate byte daysfromintubationtotrach = 1
    generate str80 attempttoweanpriortotrachyorn = "N"
    generate str80 attempttoweanoffventpriortotrans = "N"
    generate str40 didpatientdeveloppneumoniapriort = "N"
    generate byte column1 = .
    generate str8 osh = "N"
    generate str20 admitfrom = "TEST"
    generate byte daysfrominjurytoadmissiontorehab = 5
    generate str40 asiaclassificationatrehab = "C4 AIS A"
    replace asiaclassificationatrehab = "C8 AIS D" in 2
    generate str40 weanoffventduringthedayandcontin = "N"
    generate str40 ifyestopreviousquestionhowmanyda = "1"
    replace ifyestopreviousquestionhowmanyda = "PRIOR TO REHAB" in 2
    generate str40 weanoffventforfull24hoursyorn = "N"
    generate str40 v23 = "2"
    replace v23 = "1 WEANED OFF CPAP AT NIGHT" in 2
    generate str20 didpatientdecanulateyorn = "N"
    generate byte daysfromadmissiontorehabtodecanu = 2
    generate str40 didtheytransfertoicufromrehabyor = "N"
    generate str60 iftheytransferredtoicuwascausere = "N"
    generate str40 didtheydeveloppnumoniaatrehabyor = "N"
    generate byte daystodischargefromrehab = 10
    generate str40 dischargelocation = "HOME"
    generate str20 completerecord = "YES"
    generate str20 comments = "TEST"
    generate str8 deceased = "N"
    generate str16 tod = "01/02/2021"
    generate int daysfromdischargetodeath = 1
    generate str20 causeofdeath = "TEST"

    if "`variant'" == "reordered" {
        order age pat_mrn_id
    }
    else if "`variant'" == "missing_field" {
        drop causeofdeath
    }
    else if "`variant'" == "extra_field" {
        generate byte unexpected_field = 1
    }
    else if "`variant'" == "wrong_type" {
        tostring age, replace
        replace age = "NOT_NUMERIC" in 2
    }
    else if "`variant'" == "missing_required" {
        replace pat_mrn_id = . in 2
    }
    else if "`variant'" == "duplicate_id" {
        replace pat_mrn_id = pat_mrn_id[1] in 2
    }
    else if "`variant'" == "invalid_date" {
        replace dateofinjury = "NOT_A_DATE" in 2
    }
    else if "`variant'" == "negative_timing" {
        replace daystodischargefromrehab = -1 in 2
    }
    else if "`variant'" == "whitespace_blank" {
        replace sex = "   " in 2
        replace dateofinjury = "   " in 2
        replace ifyestopreviousquestionhowmanyda = "   " in 2
    }
    else if inlist("`variant'", "unknown_category", "sentinel") {
        replace sex = "NRH_SENTINEL_SECRET" in 2
    }

    export delimited using "`using'", replace
end

do "code/lib/data_contracts.do"

tempfile nrh_source nrh_log

* A conforming synthetic source passes strict mode.
nrh_test_write_fixture using "`nrh_source'", variant("valid")
nrh_validate_source_contract using "`nrh_source'", ///
    schema("validation/source_schema.csv") ///
    rules("validation/validation_rules.csv") ///
    mode("strict") log("`nrh_log'")
assert r(failed_rules) == 0
assert r(warning_rules) == 0
assert r(contract_version) == 2

* Surrounding-whitespace normalization makes whitespace-only strings blank.
nrh_test_write_fixture using "`nrh_source'", variant("whitespace_blank")
nrh_validate_source_contract using "`nrh_source'", ///
    schema("validation/source_schema.csv") ///
    rules("validation/validation_rules.csv") ///
    mode("strict") log("`nrh_log'")
assert r(failed_rules) == 0
assert r(warning_rules) == 0

* Structural failures remain fatal in both modes.
foreach nrh_variant in reordered missing_field extra_field wrong_type {
    nrh_test_write_fixture using "`nrh_source'", variant("`nrh_variant'")
    capture noisily nrh_validate_source_contract using "`nrh_source'", ///
        schema("validation/source_schema.csv") ///
        rules("validation/validation_rules.csv") ///
        mode("strict") log("`nrh_log'")
    assert _rc == 459
}

nrh_test_write_fixture using "`nrh_source'", variant("reordered")
capture noisily nrh_validate_source_contract using "`nrh_source'", ///
    schema("validation/source_schema.csv") ///
    rules("validation/validation_rules.csv") ///
    mode("development") log("`nrh_log'")
assert _rc == 459

* Content and domain failures fail strict mode.
foreach nrh_variant in missing_required duplicate_id invalid_date negative_timing unknown_category {
    nrh_test_write_fixture using "`nrh_source'", variant("`nrh_variant'")
    capture noisily nrh_validate_source_contract using "`nrh_source'", ///
        schema("validation/source_schema.csv") ///
        rules("validation/validation_rules.csv") ///
        mode("strict") log("`nrh_log'")
    assert _rc == 459
}

* Development mode warns for content drift but does not bypass structure.
nrh_test_write_fixture using "`nrh_source'", variant("unknown_category")
nrh_validate_source_contract using "`nrh_source'", ///
    schema("validation/source_schema.csv") ///
    rules("validation/validation_rules.csv") ///
    mode("development") log("`nrh_log'")
assert r(failed_rules) == 0
assert r(warning_rules) > 0

* The aggregate log must not echo the invalid source value or a row number.
nrh_test_write_fixture using "`nrh_source'", variant("sentinel")
capture noisily nrh_validate_source_contract using "`nrh_source'", ///
    schema("validation/source_schema.csv") ///
    rules("validation/validation_rules.csv") ///
    mode("strict") log("`nrh_log'")
assert _rc == 459

tempname nrh_handle
file open `nrh_handle' using "`nrh_log'", read text
file read `nrh_handle' nrh_line
while r(eof) == 0 {
    assert strpos("`nrh_line'", "NRH_SENTINEL_SECRET") == 0
    assert strpos("`nrh_line'", "row_number") == 0
    file read `nrh_handle' nrh_line
}
file close `nrh_handle'

if strtrim("`status_file'") != "" {
    tempname nrh_status
    file open `nrh_status' using "`status_file'", write text replace
    file write `nrh_status' "0" _n
    file close `nrh_status'
}

di as result "NRH source-contract synthetic tests passed."

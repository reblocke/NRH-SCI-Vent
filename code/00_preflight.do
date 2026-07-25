* Verify the project-controlled Stata dependency environment before data access.

version 18
args dependency_log

if strtrim("`dependency_log'") == "" {
    local dependency_log "dependency_preflight.log"
}

foreach required_file in ///
    "vendor/stata/manifest.csv" ///
    "vendor/stata/LICENSES.md" ///
    "vendor/stata/licenses/GPL-3.0-only.txt" ///
    "vendor/stata/licenses/MIT-coefplot.txt" ///
    "vendor/stata/licenses/MIT-schemepack.txt" {
    capture confirm file "`required_file'"
    if _rc {
        di as err "Dependency preflight could not find `required_file'."
        exit 601
    }
}

capture adopath - "vendor/stata/plus"
quietly adopath ++ "vendor/stata/plus"
discard

capture program drop _nrh_sha256_file
program define _nrh_sha256_file, rclass
    version 18
    syntax using/

    local target `"`using'"'
    tempfile sha_output
    capture erase "`sha_output'"

    if "`c(os)'" == "MacOSX" {
        capture quietly shell /usr/bin/shasum -a 256 "`target'" > "`sha_output'"
    }
    else if "`c(os)'" == "Unix" {
        capture quietly shell sha256sum "`target'" > "`sha_output'"
    }
    else if "`c(os)'" == "Windows" {
        capture quietly shell certutil -hashfile "`target'" SHA256 > "`sha_output'"
    }
    else {
        di as err "Unsupported operating system for dependency SHA-256 verification: `c(os)'."
        exit 198
    }

    capture confirm file "`sha_output'"
    if _rc {
        di as err "The operating-system SHA-256 command did not create verification output."
        exit 601
    }

    tempname sha_input
    file open `sha_input' using "`sha_output'", read text
    local actual_sha256 ""
    file read `sha_input' sha_line
    while r(eof) == 0 {
        local normalized_line = lower(strtrim(`"`sha_line'"'))
        local first_token : word 1 of `normalized_line'
        local compact_line = subinstr("`normalized_line'", " ", "", .)

        if strlen("`first_token'") == 64 & ///
            regexm("`first_token'", "^[0-9a-f][0-9a-f]*$") {
            local actual_sha256 "`first_token'"
        }
        if strlen("`compact_line'") == 64 & ///
            regexm("`compact_line'", "^[0-9a-f][0-9a-f]*$") {
            local actual_sha256 "`compact_line'"
        }
        file read `sha_input' sha_line
    }
    file close `sha_input'

    if strlen("`actual_sha256'") != 64 {
        di as err "Could not parse a 64-character SHA-256 digest."
        exit 459
    }

    return local sha256 "`actual_sha256'"
end

tempname preflight_log
capture file open `preflight_log' using "`dependency_log'", write text replace
if _rc {
    di as err "Could not create the dependency preflight log."
    capture program drop _nrh_sha256_file
    exit 603
}

file write `preflight_log' "NRH-SCI-Vent dependency preflight" _n
file write `preflight_log' "manifest=vendor/stata/manifest.csv" _n
file write `preflight_log' "vendor_path=vendor/stata/plus" _n
file write `preflight_log' "stata_version=`c(stata_version)'" _n
file write `preflight_log' "stata_edition=`c(edition_real)'" _n
file write `preflight_log' "stata_build_date=`c(born_date)'" _n
file write `preflight_log' "operating_system=`c(os)'" _n

capture quietly _nrh_sha256_file using "vendor/stata/manifest.csv"
local manifest_hash_rc = _rc
if `manifest_hash_rc' {
    file write `preflight_log' "status=failed" _n
    file write `preflight_log' "failed_check=manifest_sha256" _n
    file close `preflight_log'
    capture program drop _nrh_sha256_file
    exit `manifest_hash_rc'
}
local dependency_manifest_sha256 "`r(sha256)'"
file write `preflight_log' "manifest_sha256=`dependency_manifest_sha256'" _n

capture quietly import delimited using "vendor/stata/manifest.csv", ///
    clear varnames(1) stringcols(_all) bindquote(strict)
local manifest_import_rc = _rc
if `manifest_import_rc' {
    file write `preflight_log' "status=failed" _n
    file write `preflight_log' "failed_check=manifest_parse" _n
    file close `preflight_log'
    capture program drop _nrh_sha256_file
    exit `manifest_import_rc'
}

capture confirm variable command package version_or_date source_url license ///
    redistributable file_path sha256 verified_with_stata notes
if _rc {
    file write `preflight_log' "status=failed" _n
    file write `preflight_log' "failed_check=manifest_columns" _n
    file close `preflight_log'
    capture program drop _nrh_sha256_file
    exit 459
}

if _N != 7 {
    file write `preflight_log' "status=failed" _n
    file write `preflight_log' "failed_check=manifest_row_count" _n
    file close `preflight_log'
    capture program drop _nrh_sha256_file
    exit 459
}

sort command
capture assert command != command[_n - 1] if _n > 1
if _rc {
    file write `preflight_log' "status=failed" _n
    file write `preflight_log' "failed_check=duplicate_command" _n
    file close `preflight_log'
    capture program drop _nrh_sha256_file
    exit 459
}
sort file_path
capture assert file_path != file_path[_n - 1] if _n > 1
if _rc {
    file write `preflight_log' "status=failed" _n
    file write `preflight_log' "failed_check=duplicate_file_path" _n
    file close `preflight_log'
    capture program drop _nrh_sha256_file
    exit 459
}
capture assert redistributable == "true"
if _rc {
    file write `preflight_log' "status=failed" _n
    file write `preflight_log' "failed_check=redistribution_flag" _n
    file close `preflight_log'
    capture program drop _nrh_sha256_file
    exit 459
}

forvalues dependency_index = 1/`=_N' {
    local dependency_command = command[`dependency_index']
    local dependency_file = file_path[`dependency_index']
    local expected_sha256 = lower(sha256[`dependency_index'])

    if !regexm("`dependency_file'", ///
        "^vendor/stata/plus/[a-z]/[A-Za-z0-9_.-][A-Za-z0-9_.-]*$") | ///
        strpos("`dependency_file'", "..") {
        file write `preflight_log' "status=failed" _n
        file write `preflight_log' "failed_check=unsafe_manifest_path" _n
        file close `preflight_log'
        capture program drop _nrh_sha256_file
        exit 198
    }
    if strlen("`expected_sha256'") != 64 | ///
        !regexm("`expected_sha256'", "^[0-9a-f][0-9a-f]*$") {
        file write `preflight_log' "status=failed" _n
        file write `preflight_log' "failed_check=invalid_expected_sha256" _n
        file close `preflight_log'
        capture program drop _nrh_sha256_file
        exit 459
    }

    capture confirm file "`dependency_file'"
    if _rc {
        file write `preflight_log' "status=failed" _n
        file write `preflight_log' "failed_check=missing_dependency_file" _n
        file write `preflight_log' "file=`dependency_file'" _n
        file close `preflight_log'
        capture program drop _nrh_sha256_file
        exit 601
    }

    capture quietly _nrh_sha256_file using "`dependency_file'"
    local dependency_hash_rc = _rc
    if `dependency_hash_rc' {
        file write `preflight_log' "status=failed" _n
        file write `preflight_log' "failed_check=dependency_sha256" _n
        file write `preflight_log' "file=`dependency_file'" _n
        file close `preflight_log'
        capture program drop _nrh_sha256_file
        exit `dependency_hash_rc'
    }
    local actual_sha256 "`r(sha256)'"
    if "`actual_sha256'" != "`expected_sha256'" {
        file write `preflight_log' "status=failed" _n
        file write `preflight_log' "failed_check=dependency_hash_mismatch" _n
        file write `preflight_log' "file=`dependency_file'" _n
        file write `preflight_log' "expected_sha256=`expected_sha256'" _n
        file write `preflight_log' "actual_sha256=`actual_sha256'" _n
        file close `preflight_log'
        capture program drop _nrh_sha256_file
        exit 459
    }

    if substr("`dependency_command'", 1, 7) == "scheme:" {
        local scheme_name = substr("`dependency_command'", 8, .)
        capture quietly findfile "scheme-`scheme_name'.scheme"
    }
    else {
        capture quietly which `dependency_command'
        if !_rc {
            capture quietly findfile "`dependency_command'.ado"
        }
    }
    local dependency_resolve_rc = _rc
    if `dependency_resolve_rc' {
        file write `preflight_log' "status=failed" _n
        file write `preflight_log' "failed_check=dependency_resolution" _n
        file write `preflight_log' "dependency=`dependency_command'" _n
        file close `preflight_log'
        capture program drop _nrh_sha256_file
        exit `dependency_resolve_rc'
    }

    local resolved_file = lower(subinstr(strtrim(`"`r(fn)'"'), "\", "/", .))
    local expected_file = lower(subinstr("`c(pwd)'/`dependency_file'", "\", "/", .))
    if substr("`resolved_file'", 1, 1) != "/" & ///
        !regexm("`resolved_file'", "^[a-z]:/") {
        local resolved_file = lower(subinstr("`c(pwd)'/`resolved_file'", "\", "/", .))
    }
    while strpos("`resolved_file'", "//") {
        local resolved_file = subinstr("`resolved_file'", "//", "/", .)
    }
    while strpos("`expected_file'", "//") {
        local expected_file = subinstr("`expected_file'", "//", "/", .)
    }

    if "`resolved_file'" != "`expected_file'" {
        file write `preflight_log' "status=failed" _n
        file write `preflight_log' "failed_check=uncontrolled_dependency_resolution" _n
        file write `preflight_log' "dependency=`dependency_command'" _n
        file close `preflight_log'
        capture program drop _nrh_sha256_file
        exit 459
    }

    file write `preflight_log' "dependency=`dependency_command'" _n
    file write `preflight_log' "resolved_file=`dependency_file'" _n
    file write `preflight_log' "expected_sha256=`expected_sha256'" _n
    file write `preflight_log' "actual_sha256=`actual_sha256'" _n
}

file write `preflight_log' "status=passed" _n
file close `preflight_log'

global NRH_DEPENDENCY_MANIFEST_SHA256 "`dependency_manifest_sha256'"
capture program drop _nrh_sha256_file

di as txt "NRH dependency preflight passed."

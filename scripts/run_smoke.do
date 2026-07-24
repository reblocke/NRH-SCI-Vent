* Batch adapter for smoke launchers. Stata for macOS may return process status
* zero even when a do-file returns an error, so write the Stata return code to
* a launcher-owned temporary sidecar.

version 18
clear
set more off

args repo_root status_file data_dir output_root run_id

local nrh_wrapper_rc 0
capture confirm file "`repo_root'/run_all.do"
if _rc {
    di as err "Could not find run_all.do under the launcher repository root."
    local nrh_wrapper_rc 601
}
else {
    quietly cd "`repo_root'"
    capture noisily do run_all.do smoke ///
        "`data_dir'" "`output_root'" "`run_id'"
    local nrh_wrapper_rc = _rc
}

tempname nrh_status
capture file open `nrh_status' using "`status_file'", write text replace
if _rc {
    di as err "Could not write the launcher status sidecar."
    exit, clear STATA
}
file write `nrh_status' "`nrh_wrapper_rc'" _n
file close `nrh_status'

exit, clear STATA

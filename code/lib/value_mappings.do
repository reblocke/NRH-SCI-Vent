* Contract-backed categorical value mappings for NRH-SCI-Vent.
*
* The mapper normalizes source strings with Unicode NFC, surrounding-space
* trimming, and uppercase conversion. It emits aggregate diagnostics only:
* target variable, mapped count, missing count, rejected count, and status.
* It never emits source literals, row numbers, identifiers, or source paths.

version 18

capture program drop nrh_apply_value_mapping
program define nrh_apply_value_mapping, rclass
    version 18
    syntax varname(string), SOURCEFIELD(string) TARGET(name) ///
        CONTRACT(string) DIAGNOSTICS(string)

    local nrh_source_field = strtrim("`sourcefield'")
    local nrh_target = strtrim("`target'")
    if "`nrh_source_field'" == "" | "`nrh_target'" == "" {
        di as err "The value-mapping source field and target variable are required."
        exit 198
    }

    capture confirm file "`contract'"
    if _rc {
        di as err "The public value-mapping contract is unavailable."
        exit 601
    }

    capture confirm new variable `nrh_target'
    if _rc {
        di as err "The requested value-mapping target already exists."
        exit 110
    }

    preserve
    capture quietly import delimited using "`contract'", clear varnames(1) ///
        stringcols(_all) bindquote(strict) encoding("UTF-8")
    if _rc {
        restore
        di as err "The public value-mapping contract could not be parsed."
        exit 459
    }

    capture confirm variable mapping_id source_field input_kind ///
        normalized_literal target_variable action target_code target_label ///
        decision_id unlisted_action notes
    if _rc {
        restore
        di as err "The public value-mapping contract is malformed."
        exit 459
    }

    quietly keep if source_field == "`nrh_source_field'" & ///
        target_variable == "`nrh_target'"
    quietly count
    local nrh_mapping_count = r(N)
    if `nrh_mapping_count' == 0 {
        restore
        di as err "No approved mapping exists for the requested source and target."
        exit 459
    }

    quietly count if strtrim(mapping_id) == "" | ///
        strtrim(decision_id) == "" | ///
        !inlist(input_kind, "literal", "blank") | ///
        !inlist(action, "assign", "set_missing") | ///
        unlisted_action != "reject"
    if r(N) {
        restore
        di as err "The requested value-mapping rows have invalid metadata."
        exit 459
    }

    quietly count if ///
        (input_kind == "blank" & ///
            (normalized_literal != "" | action != "set_missing")) | ///
        (input_kind == "literal" & normalized_literal == "") | ///
        (input_kind == "literal" & normalized_literal != ///
            ustrupper(ustrtrim(ustrnormalize(normalized_literal, "nfc"))))
    if r(N) {
        restore
        di as err "The requested value-mapping rows violate normalization rules."
        exit 459
    }

    tempvar nrh_code
    quietly generate double `nrh_code' = real(target_code)
    quietly count if ///
        (action == "assign" & ///
            (missing(`nrh_code') | `nrh_code' != floor(`nrh_code') | ///
                `nrh_code' < 0 | strtrim(target_label) == "")) | ///
        (action == "set_missing" & ///
            (strtrim(target_code) != "" | strtrim(target_label) != ""))
    if r(N) {
        restore
        di as err "The requested value-mapping rows have invalid codes or labels."
        exit 459
    }

    tempvar nrh_duplicate_id
    quietly duplicates tag mapping_id, generate(`nrh_duplicate_id')
    quietly count if `nrh_duplicate_id' > 0
    if r(N) {
        restore
        di as err "The requested value-mapping rows contain duplicate mapping IDs."
        exit 459
    }

    tempvar nrh_duplicate
    quietly duplicates tag input_kind normalized_literal, ///
        generate(`nrh_duplicate')
    quietly count if `nrh_duplicate' > 0
    if r(N) {
        restore
        di as err "The requested value-mapping rows contain conflicting inputs."
        exit 459
    }

    quietly count if input_kind == "blank"
    if r(N) != 1 {
        restore
        di as err "The requested value mapping must define exactly one blank action."
        exit 459
    }

    * A target code must always have one stable value label.
    tempvar nrh_label_conflict
    quietly bysort `nrh_code' (target_label): generate byte ///
        `nrh_label_conflict' = ///
        action == "assign" & target_label[1] != target_label[_N]
    quietly count if `nrh_label_conflict'
    if r(N) {
        restore
        di as err "The requested value mapping has conflicting target labels."
        exit 459
    }

    forvalues nrh_i = 1/`nrh_mapping_count' {
        local nrh_kind_`nrh_i' `"`=input_kind[`nrh_i']'"'
        local nrh_literal_`nrh_i' `"`=normalized_literal[`nrh_i']'"'
        local nrh_action_`nrh_i' `"`=action[`nrh_i']'"'
        local nrh_code_`nrh_i' = `nrh_code'[`nrh_i']
        local nrh_label_`nrh_i' `"`=target_label[`nrh_i']'"'
    }
    restore

    tempvar nrh_normalized nrh_matched
    quietly generate strL `nrh_normalized' = ///
        ustrupper(ustrtrim(ustrnormalize(`varlist', "nfc")))
    quietly generate byte `nrh_matched' = 0
    quietly generate long `nrh_target' = .

    local nrh_label_codes ""
    local nrh_label_count 0
    forvalues nrh_i = 1/`nrh_mapping_count' {
        local nrh_kind `"`nrh_kind_`nrh_i''"'
        local nrh_literal `"`nrh_literal_`nrh_i''"'
        local nrh_action "`nrh_action_`nrh_i''"

        if "`nrh_kind'" == "blank" {
            quietly replace `nrh_matched' = 1 if missing(`nrh_normalized')
        }
        else {
            if "`nrh_action'" == "assign" {
                quietly replace `nrh_target' = `nrh_code_`nrh_i'' if ///
                    `nrh_normalized' == `"`nrh_literal'"'
            }
            quietly replace `nrh_matched' = 1 if ///
                `nrh_normalized' == `"`nrh_literal'"'
        }

        if "`nrh_action'" == "assign" {
            local nrh_code = `nrh_code_`nrh_i''
            if strpos(" `nrh_label_codes' ", " `nrh_code' ") == 0 {
                local nrh_label_count = `nrh_label_count' + 1
                local nrh_label_code_`nrh_label_count' = `nrh_code'
                local nrh_label_text_`nrh_label_count' ///
                    `"`nrh_label_`nrh_i''"'
                local nrh_label_codes "`nrh_label_codes' `nrh_code'"
            }
        }
    }

    quietly count if `nrh_matched' & !missing(`nrh_target')
    local nrh_mapped_count = r(N)
    quietly count if `nrh_matched' & missing(`nrh_target')
    local nrh_missing_count = r(N)
    quietly count if !`nrh_matched'
    local nrh_rejected_count = r(N)
    local nrh_status = cond(`nrh_rejected_count', "fail", "pass")

    tempname nrh_diagnostics
    capture confirm file "`diagnostics'"
    if _rc {
        capture file open `nrh_diagnostics' using "`diagnostics'", ///
            write text replace
        if _rc {
            quietly drop `nrh_target'
            di as err "The aggregate value-mapping diagnostics could not be created."
            exit 603
        }
        file write `nrh_diagnostics' ///
            "target_variable,mapped_count,missing_count,rejected_count,status" _n
    }
    else {
        capture file open `nrh_diagnostics' using "`diagnostics'", ///
            write text append
        if _rc {
            quietly drop `nrh_target'
            di as err "The aggregate value-mapping diagnostics could not be updated."
            exit 603
        }
    }
    file write `nrh_diagnostics' ///
        "`nrh_target',`nrh_mapped_count',`nrh_missing_count'," ///
        "`nrh_rejected_count',`nrh_status'" _n
    file close `nrh_diagnostics'

    if `nrh_rejected_count' {
        quietly drop `nrh_target'
        di as err "An unlisted value was rejected by the approved mapping contract."
        exit 459
    }

    local nrh_value_label "nrh_`nrh_target'_label"
    forvalues nrh_i = 1/`nrh_label_count' {
        if `nrh_i' == 1 {
            label define `nrh_value_label' ///
                `nrh_label_code_`nrh_i'' ///
                `"`nrh_label_text_`nrh_i''"', replace
        }
        else {
            label define `nrh_value_label' ///
                `nrh_label_code_`nrh_i'' ///
                `"`nrh_label_text_`nrh_i''"', add
        }
    }
    label values `nrh_target' `nrh_value_label'

    return scalar mapped_count = `nrh_mapped_count'
    return scalar missing_count = `nrh_missing_count'
    return scalar rejected_count = `nrh_rejected_count'
    return local target_variable "`nrh_target'"
end

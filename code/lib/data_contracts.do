* Public source-contract version 2 validator for NRH-SCI-Vent.
*
* The validator emits only aggregate invalid counts. It never lists source
* values, identifiers, row numbers, or input paths.

version 18

capture program drop nrh_validate_source_contract
program define nrh_validate_source_contract, rclass
    version 18
    syntax using/, SCHEMA(string) RULES(string) MODE(string) LOG(string)

    local nrh_mode = lower(strtrim("`mode'"))
    if !inlist("`nrh_mode'", "strict", "development") {
        di as err "Source-contract mode must be strict or development."
        exit 198
    }

    foreach nrh_contract_file in "`schema'" "`rules'" {
        capture confirm file "`nrh_contract_file'"
        if _rc {
            di as err "A required public source-contract file is unavailable."
            exit 601
        }
    }

    capture confirm file "`using'"
    if _rc {
        di as err "The required source CSV is unavailable."
        exit 601
    }

    quietly import delimited using "`schema'", clear varnames(1) ///
        stringcols(_all) bindquote(strict) encoding("UTF-8")
    capture confirm variable source_field required storage_type semantic_type ///
        allowed_values missing_allowed sensitivity used_by notes
    if _rc {
        di as err "The public source-schema contract is malformed."
        exit 459
    }

    quietly count
    local nrh_field_count = r(N)
    if `nrh_field_count' == 0 {
        di as err "The public source-schema contract has no fields."
        exit 459
    }

    quietly count if strtrim(source_field) == ""
    if r(N) {
        di as err "The public source-schema contract has a blank field name."
        exit 459
    }
    quietly duplicates tag source_field, generate(__nrh_duplicate)
    quietly count if __nrh_duplicate
    if r(N) {
        di as err "The public source-schema contract has duplicate field names."
        exit 459
    }
    drop __nrh_duplicate

    quietly count if !inlist(required, "true", "false") | ///
        !inlist(missing_allowed, "true", "false")
    if r(N) {
        di as err "The public source-schema contract has an invalid boolean."
        exit 459
    }
    quietly count if required != "true"
    if r(N) {
        di as err "Every field in the ordered source contract must be required."
        exit 459
    }
    quietly count if !inlist(storage_type, "numeric", "string")
    if r(N) {
        di as err "The public source-schema contract has an invalid storage type."
        exit 459
    }
    quietly count if !inlist(semantic_type, "identifier", "integer", ///
        "date_mdy", "categorical", "numeric_with_exceptions", ///
        "free_text", "empty_placeholder")
    if r(N) {
        di as err "The public source-schema contract has an invalid semantic type."
        exit 459
    }
    quietly count if !inlist(sensitivity, "direct_identifier", ///
        "quasi_identifier", "sensitive_clinical", ///
        "restricted_free_text", "operational")
    if r(N) {
        di as err "The public source-schema contract has an invalid sensitivity class."
        exit 459
    }
    quietly count if semantic_type == "categorical" & strtrim(allowed_values) == ""
    if r(N) {
        di as err "A categorical source field lacks an allowed-value domain."
        exit 459
    }
    quietly count if semantic_type == "free_text" & strtrim(allowed_values) != ""
    if r(N) {
        di as err "A free-text source field may not publish an allowed-value domain."
        exit 459
    }
    quietly count if ///
        (inlist(semantic_type, "identifier", "integer", "empty_placeholder") & ///
            storage_type != "numeric") | ///
        (inlist(semantic_type, "date_mdy", "categorical", ///
            "numeric_with_exceptions", "free_text") & storage_type != "string")
    if r(N) {
        di as err "A source field has incompatible storage and semantic types."
        exit 459
    }

    forvalues nrh_i = 1/`nrh_field_count' {
        local nrh_field_`nrh_i' = source_field[`nrh_i']
        local nrh_required_`nrh_i' = required[`nrh_i']
        local nrh_storage_`nrh_i' = storage_type[`nrh_i']
        local nrh_semantic_`nrh_i' = semantic_type[`nrh_i']
        local nrh_allowed_`nrh_i' = allowed_values[`nrh_i']
        local nrh_missing_`nrh_i' = missing_allowed[`nrh_i']
    }

    quietly import delimited using "`rules'", clear varnames(1) ///
        stringcols(_all) bindquote(strict) encoding("UTF-8")
    capture confirm variable rule_id source_field rule_type parameters ///
        strict_action development_action disclosure_output notes
    if _rc {
        di as err "The public validation-rules contract is malformed."
        exit 459
    }

    quietly count
    local nrh_rule_count = r(N)
    if `nrh_rule_count' == 0 {
        di as err "The public validation-rules contract has no rules."
        exit 459
    }
    quietly count if strtrim(rule_id) == "" | strtrim(source_field) == ""
    if r(N) {
        di as err "The public validation-rules contract has a blank identifier."
        exit 459
    }
    quietly duplicates tag rule_id, generate(__nrh_duplicate)
    quietly count if __nrh_duplicate
    if r(N) {
        di as err "The public validation-rules contract has duplicate rule IDs."
        exit 459
    }
    drop __nrh_duplicate
    quietly count if !inlist(rule_type, "unique_nonmissing", ///
        "parseable_date", "integer", "nonnegative", "all_missing", ///
        "numeric_or_allowed_literal")
    if r(N) {
        di as err "The public validation-rules contract has an invalid rule type."
        exit 459
    }
    quietly count if strict_action != "fail" | ///
        !inlist(development_action, "fail", "warn") | ///
        disclosure_output != "aggregate_count_only"
    if r(N) {
        di as err "The public validation-rules contract has an invalid action or disclosure policy."
        exit 459
    }
    quietly count if ///
        (rule_type == "parseable_date" & parameters != "mask=MD20Y") | ///
        (rule_type == "numeric_or_allowed_literal" & ///
            parameters != "nonnegative=true") | ///
        (!inlist(rule_type, "parseable_date", ///
            "numeric_or_allowed_literal") & strtrim(parameters) != "")
    if r(N) {
        di as err "The public validation-rules contract has invalid parameters."
        exit 459
    }

    generate byte __nrh_known_field = 0
    forvalues nrh_i = 1/`nrh_field_count' {
        quietly replace __nrh_known_field = 1 if ///
            source_field == "`nrh_field_`nrh_i''"
    }
    quietly count if __nrh_known_field == 0
    if r(N) {
        di as err "A validation rule references an unknown source field."
        exit 459
    }
    drop __nrh_known_field

    forvalues nrh_i = 1/`nrh_rule_count' {
        local nrh_rule_id_`nrh_i' = rule_id[`nrh_i']
        local nrh_rule_field_`nrh_i' = source_field[`nrh_i']
        local nrh_rule_type_`nrh_i' = rule_type[`nrh_i']
        local nrh_rule_parameters_`nrh_i' = parameters[`nrh_i']
        local nrh_rule_dev_`nrh_i' = development_action[`nrh_i']
    }

    tempname nrh_contract_log
    capture file open `nrh_contract_log' using "`log'", write text replace
    if _rc {
        di as err "The aggregate source-contract log could not be created."
        exit 603
    }
    file write `nrh_contract_log' ///
        "rule_id,source_field,status,invalid_count" _n

    capture quietly import delimited using "`using'", clear varnames(1) ///
        case(lower) bindquote(strict) encoding("UTF-8")
    if _rc {
        file write `nrh_contract_log' "schema.csv_parse,_source,fail,1" _n
        file close `nrh_contract_log'
        di as err "The source CSV could not be parsed."
        exit 459
    }

    local nrh_structural_failures 0
    local nrh_content_failures 0
    local nrh_warnings 0

    local nrh_actual_count = c(k)
    local nrh_count_invalid = (`nrh_actual_count' != `nrh_field_count')
    local nrh_count_status = cond(`nrh_count_invalid', "fail", "pass")
    file write `nrh_contract_log' ///
        "schema.column_count,_source,`nrh_count_status',`nrh_count_invalid'" _n
    if `nrh_count_invalid' local nrh_structural_failures = `nrh_structural_failures' + 1

    unab nrh_actual_fields : _all
    local nrh_expected_fields ""
    forvalues nrh_i = 1/`nrh_field_count' {
        local nrh_expected_fields ///
            "`nrh_expected_fields' `nrh_field_`nrh_i''"
    }
    local nrh_expected_fields = strtrim("`nrh_expected_fields'")
    local nrh_order_invalid = ("`nrh_actual_fields'" != "`nrh_expected_fields'")
    local nrh_order_status = cond(`nrh_order_invalid', "fail", "pass")
    file write `nrh_contract_log' ///
        "schema.field_order,_source,`nrh_order_status',`nrh_order_invalid'" _n
    if `nrh_order_invalid' local nrh_structural_failures = `nrh_structural_failures' + 1

    if `nrh_count_invalid' == 0 & `nrh_order_invalid' == 0 {
        forvalues nrh_i = 1/`nrh_field_count' {
            local nrh_field "`nrh_field_`nrh_i''"
            local nrh_storage "`nrh_storage_`nrh_i''"
            local nrh_type_invalid 0
            if "`nrh_storage'" == "numeric" {
                capture confirm numeric variable `nrh_field'
                if _rc local nrh_type_invalid 1
            }
            else {
                capture confirm string variable `nrh_field'
                if _rc local nrh_type_invalid 1
            }
            local nrh_type_status = cond(`nrh_type_invalid', "fail", "pass")
            file write `nrh_contract_log' ///
                "schema.`nrh_field'.storage,`nrh_field',`nrh_type_status',`nrh_type_invalid'" _n
            if `nrh_type_invalid' local nrh_structural_failures = `nrh_structural_failures' + 1
        }
    }

    * Structural failures remain fatal in development mode.
    if `nrh_structural_failures' {
        file close `nrh_contract_log'
        return scalar failed_rules = `nrh_structural_failures'
        return scalar warning_rules = 0
        return scalar contract_version = 2
        di as err "The source CSV failed structural contract validation."
        exit 459
    }

    * Schema-driven content checks.
    forvalues nrh_i = 1/`nrh_field_count' {
        local nrh_field "`nrh_field_`nrh_i''"
        local nrh_storage "`nrh_storage_`nrh_i''"
        local nrh_semantic "`nrh_semantic_`nrh_i''"
        local nrh_allowed `"`nrh_allowed_`nrh_i''"'
        local nrh_missing "`nrh_missing_`nrh_i''"

        if "`nrh_missing'" == "false" {
            quietly count if missing(`nrh_field')
            local nrh_invalid = r(N)
            local nrh_status = cond(`nrh_invalid', ///
                cond("`nrh_mode'" == "strict", "fail", "warn"), "pass")
            file write `nrh_contract_log' ///
                "schema.`nrh_field'.missing,`nrh_field',`nrh_status',`nrh_invalid'" _n
            if `nrh_invalid' {
                if "`nrh_mode'" == "strict" local nrh_content_failures = `nrh_content_failures' + 1
                else local nrh_warnings = `nrh_warnings' + 1
            }
        }

        if "`nrh_semantic'" == "integer" {
            quietly count if !missing(`nrh_field') & `nrh_field' != floor(`nrh_field')
            local nrh_invalid = r(N)
            local nrh_status = cond(`nrh_invalid', ///
                cond("`nrh_mode'" == "strict", "fail", "warn"), "pass")
            file write `nrh_contract_log' ///
                "schema.`nrh_field'.integer,`nrh_field',`nrh_status',`nrh_invalid'" _n
            if `nrh_invalid' {
                if "`nrh_mode'" == "strict" local nrh_content_failures = `nrh_content_failures' + 1
                else local nrh_warnings = `nrh_warnings' + 1
            }
        }

        if "`nrh_semantic'" == "categorical" {
            tempvar nrh_normalized
            quietly generate strL `nrh_normalized' = ///
                ustrupper(ustrtrim(ustrnormalize(`nrh_field', "nfc")))
            quietly count if !missing(`nrh_normalized') & ///
                strpos(" | `nrh_allowed' | ", " | " + `nrh_normalized' + " | ") == 0
            local nrh_invalid = r(N)
            local nrh_status = cond(`nrh_invalid', ///
                cond("`nrh_mode'" == "strict", "fail", "warn"), "pass")
            file write `nrh_contract_log' ///
                "schema.`nrh_field'.allowed,`nrh_field',`nrh_status',`nrh_invalid'" _n
            if `nrh_invalid' {
                if "`nrh_mode'" == "strict" local nrh_content_failures = `nrh_content_failures' + 1
                else local nrh_warnings = `nrh_warnings' + 1
            }
        }
    }

    * Additional rule checks.
    forvalues nrh_i = 1/`nrh_rule_count' {
        local nrh_rule_id "`nrh_rule_id_`nrh_i''"
        local nrh_field "`nrh_rule_field_`nrh_i''"
        local nrh_rule_type "`nrh_rule_type_`nrh_i''"
        local nrh_parameters "`nrh_rule_parameters_`nrh_i''"
        local nrh_dev_action "`nrh_rule_dev_`nrh_i''"
        local nrh_invalid 0

        if "`nrh_rule_type'" == "unique_nonmissing" {
            tempvar nrh_duplicate
            quietly duplicates tag `nrh_field' if !missing(`nrh_field'), ///
                generate(`nrh_duplicate')
            quietly count if missing(`nrh_field') | `nrh_duplicate' > 0
            local nrh_invalid = r(N)
        }
        else if "`nrh_rule_type'" == "parseable_date" {
            quietly count if !missing(ustrtrim(`nrh_field')) & ///
                date(ustrtrim(`nrh_field'), "MD20Y") == .
            local nrh_invalid = r(N)
        }
        else if "`nrh_rule_type'" == "integer" {
            quietly count if !missing(`nrh_field') & ///
                `nrh_field' != floor(`nrh_field')
            local nrh_invalid = r(N)
        }
        else if "`nrh_rule_type'" == "nonnegative" {
            quietly count if !missing(`nrh_field') & `nrh_field' < 0
            local nrh_invalid = r(N)
        }
        else if "`nrh_rule_type'" == "all_missing" {
            quietly count if !missing(`nrh_field')
            local nrh_invalid = r(N)
        }
        else if "`nrh_rule_type'" == "numeric_or_allowed_literal" {
            tempvar nrh_normalized nrh_numeric
            quietly generate strL `nrh_normalized' = ///
                ustrupper(ustrtrim(ustrnormalize(`nrh_field', "nfc")))
            quietly generate double `nrh_numeric' = real(`nrh_normalized')

            local nrh_allowed ""
            forvalues nrh_j = 1/`nrh_field_count' {
                if "`nrh_field'" == "`nrh_field_`nrh_j''" {
                    local nrh_allowed `"`nrh_allowed_`nrh_j''"'
                }
            }
            quietly count if !missing(`nrh_normalized') & ///
                ((missing(`nrh_numeric') & ///
                strpos(" | `nrh_allowed' | ", " | " + `nrh_normalized' + " | ") == 0) | ///
                (!missing(`nrh_numeric') & ///
                (`nrh_numeric' < 0 | `nrh_numeric' != floor(`nrh_numeric'))))
            local nrh_invalid = r(N)
        }

        local nrh_action "fail"
        if "`nrh_mode'" == "development" local nrh_action "`nrh_dev_action'"
        local nrh_status = cond(`nrh_invalid', "`nrh_action'", "pass")
        file write `nrh_contract_log' ///
            "`nrh_rule_id',`nrh_field',`nrh_status',`nrh_invalid'" _n
        if `nrh_invalid' {
            if "`nrh_action'" == "fail" local nrh_content_failures = `nrh_content_failures' + 1
            else local nrh_warnings = `nrh_warnings' + 1
        }
    }

    file close `nrh_contract_log'
    return scalar failed_rules = `nrh_content_failures'
    return scalar warning_rules = `nrh_warnings'
    return scalar contract_version = 2

    if `nrh_content_failures' {
        di as err "The source CSV failed content contract validation."
        exit 459
    }
end

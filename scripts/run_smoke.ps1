[Console]::Error.WriteLine(
    "The public smoke profile is retired because no synthetic cohort is authorized."
)
[Console]::Error.WriteLine(
    "Run .\scripts\run_verify.ps1 for explicitly no-data public verification."
)
exit 198

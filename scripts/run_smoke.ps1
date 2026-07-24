$ErrorActionPreference = "Stop"

$repoRoot = Split-Path -Parent $PSScriptRoot

$stataBin = if ($env:STATA_BIN) { $env:STATA_BIN } else { "stata-mp" }
$batchFlag = if ($env:STATA_BATCH_FLAG) { $env:STATA_BATCH_FLAG } else { "/e" }

if (-not (Get-Command $stataBin -ErrorAction SilentlyContinue)) {
    [Console]::Error.WriteLine("Could not find or execute Stata: $stataBin")
    exit 127
}

$tempDir = Join-Path ([System.IO.Path]::GetTempPath()) ("nrh-smoke-" + [guid]::NewGuid().ToString("N"))
$statusFile = Join-Path $tempDir "stata-status.txt"
$previousLocation = Get-Location
$exitCode = 1

New-Item -ItemType Directory -Path $tempDir | Out-Null

try {
    Set-Location $tempDir
    $stataCommand = 'do "{0}" "{1}" "{2}"' -f (Join-Path $repoRoot "scripts/run_smoke.do"), $repoRoot, $statusFile
    foreach ($stataArg in $args) {
        if ($stataArg.Contains('"')) {
            throw "Smoke-launcher arguments may not contain a double quote."
        }
        $stataCommand += ' "' + $stataArg + '"'
    }

    & $stataBin $batchFlag $stataCommand
    $stataProcessRc = $LASTEXITCODE

    if (-not (Test-Path $statusFile)) {
        [Console]::Error.WriteLine("Stata did not produce the smoke-run status sidecar.")
        if ($stataProcessRc -ne 0) {
            $exitCode = $stataProcessRc
        }
    }
    else {
        $statusText = (Get-Content -Raw $statusFile).Trim()
        $stataRunRc = 0
        if (-not [int]::TryParse($statusText, [ref]$stataRunRc)) {
            [Console]::Error.WriteLine("Stata produced an invalid smoke-run status: $statusText")
        }
        elseif ($stataRunRc -ne 0) {
            [Console]::Error.WriteLine("NRH smoke orchestration failed with Stata return code $stataRunRc.")
        }
        elseif ($stataProcessRc -ne 0) {
            [Console]::Error.WriteLine("Stata exited unexpectedly with process status $stataProcessRc.")
            $exitCode = $stataProcessRc
        }
        else {
            Write-Output "NRH smoke orchestration completed successfully."
            $exitCode = 0
        }
    }
}
finally {
    Set-Location $previousLocation
    Remove-Item -Force -ErrorAction SilentlyContinue $statusFile
    Remove-Item -Force -ErrorAction SilentlyContinue (Join-Path $tempDir "run_smoke.log")
    Remove-Item -Force -ErrorAction SilentlyContinue (Join-Path $tempDir "run_smoke.smcl")
    Remove-Item -Force -ErrorAction SilentlyContinue $tempDir
}

exit $exitCode

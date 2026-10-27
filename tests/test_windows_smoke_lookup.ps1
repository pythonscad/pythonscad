$ErrorActionPreference = 'Stop'
. (Join-Path $PSScriptRoot '../scripts/release-smoke/common-smoke.ps1')

function Assert-Equal($Actual, $Expected) {
    if ($Actual -cne $Expected) {
        throw "Expected '$Expected', got '$Actual'"
    }
}

function Assert-Throws([scriptblock]$Action, [string]$Expected) {
    try { & $Action } catch {
        if ($_.Exception.Message -notlike "*$Expected*") { throw }
        return
    }
    throw "Expected an error containing '$Expected'"
}

function gh {
    $script:ghArguments = @($args)
    $global:LASTEXITCODE = $script:ghExitCode
    $script:ghOutput
}

$script:ghExitCode = 0
$script:ghOutput = '[{"databaseId":12345}]'
foreach ($ref in @('master', 'feature/with-dashes', '')) {
    Assert-Equal (Get-LatestSuccessfulRun -Repo owner/repo -Workflow build.yml -Ref $ref) '12345'
    $arguments = $script:ghArguments
    if ($arguments -contains '--jq') { throw 'Lookup must not send jq expressions through native quoting' }
    Assert-Equal $arguments[([array]::IndexOf($arguments, '--limit') + 1)] '1'
    if ($ref) {
        Assert-Equal $arguments[([array]::IndexOf($arguments, '--branch') + 1)] $ref
    } elseif ($arguments -contains '--branch') {
        throw 'Unfiltered lookup must omit --branch'
    }
}

$script:ghOutput = '[]'
Assert-Throws { Get-LatestSuccessfulRun -Repo owner/repo -Workflow build.yml -Ref master } 'No successful build.yml run'
$script:ghOutput = ''
$script:ghExitCode = 1
Assert-Throws { Get-LatestSuccessfulRun -Repo owner/repo -Workflow build.yml -Ref master } 'Failed to query build.yml runs'
Write-Output 'PASS: Windows smoke run lookup'

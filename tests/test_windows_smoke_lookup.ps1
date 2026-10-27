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
    if ($args[0] -eq 'run' -and $args[1] -eq 'view') {
        $global:LASTEXITCODE = 0
        '{"headSha":"' + $script:actualCommit + '"}'
        return
    }
    $script:ghArguments = @($args)
    $global:LASTEXITCODE = $script:ghExitCode
    $script:ghOutput
}

$commit = 'a' * 40
$script:actualCommit = $commit
$script:ghExitCode = 0
$script:ghOutput = '[{"databaseId":12345}]'
foreach ($ref in @('master', 'feature/with-dashes', '')) {
    Assert-Equal (Get-LatestSuccessfulRun -Repo owner/repo -Workflow build.yml -Ref $ref -Commit $commit) '12345'
    $arguments = $script:ghArguments
    if ($arguments -contains '--jq') { throw 'Lookup must not send jq expressions through native quoting' }
    Assert-Equal $arguments[([array]::IndexOf($arguments, '--limit') + 1)] '1'
    Assert-Equal $arguments[([array]::IndexOf($arguments, '--commit') + 1)] $commit
    if ($ref) {
        Assert-Equal $arguments[([array]::IndexOf($arguments, '--branch') + 1)] $ref
    } elseif ($arguments -contains '--branch') {
        throw 'Unfiltered lookup must omit --branch'
    }
}

$script:ghOutput = '[]'
Assert-Throws { Get-LatestSuccessfulRun -Repo owner/repo -Workflow build.yml -Ref master -Commit $commit } 'No successful build.yml run'
$script:ghOutput = ''
$script:ghExitCode = 1
Assert-Throws { Get-LatestSuccessfulRun -Repo owner/repo -Workflow build.yml -Ref master -Commit $commit } 'Failed to query build.yml runs'
$script:ghExitCode = 0
$script:ghOutput = '[{"databaseId":12345}]'
$script:actualCommit = 'b' * 40
Assert-Throws { Get-LatestSuccessfulRun -Repo owner/repo -Workflow build.yml -Commit $commit } 'expected'
Assert-Throws { Assert-SmokeRunCommit -Repo owner/repo -RunId 12345 -Commit $commit } 'expected'
Assert-Throws { Get-ExpectedSmokeCommit 'master' } '40-character'
Write-Output 'PASS: Windows smoke run lookup'

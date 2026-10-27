param([Parameter(Mandatory = $true)][string]$PythonExecutable)

$ErrorActionPreference = 'Stop'
. (Join-Path $PSScriptRoot '../scripts/release-smoke/common-smoke.ps1')

$directory = Join-Path ([System.IO.Path]::GetTempPath()) ([System.IO.Path]::GetRandomFileName())
New-Item -ItemType Directory -Path $directory | Out-Null
$previousEncoding = [Console]::InputEncoding
try {
    # Exercise the BOM-producing console encoding used by Windows PowerShell.
    # PowerShell 7 also checks the UTF-8 source round trip on other hosts.
    [Console]::InputEncoding = [System.Text.UTF8Encoding]::new($true)
    $inputFile = Join-Path $directory 'input.py'
    $logFile = Join-Path $directory 'output.log'
    $source = "value = '" + [char]955 + "'; assert ord(value) == 955; print('UTF8_PIPE_OK')"
    Write-Utf8Content -Path $inputFile -Value $source
    $probe = 'import sys; data = sys.stdin.buffer.read(); assert not data.startswith(bytes([239,187,191])), repr(data); source = data.decode("utf-8"); assert chr(955) in source; exec(compile(source, "<stdin>", "exec"))'
    Invoke-PythonSCAD -ExecutablePath $PythonExecutable -Arguments @('-c', $probe) `
        -InputFile $inputFile -LogFile $logFile
    if ((Get-Content -LiteralPath $logFile -Raw) -notmatch 'UTF8_PIPE_OK') {
        throw 'Piped Python source did not execute'
    }
    Write-Output 'PASS: UTF-8 smoke input has no BOM and preserves non-ASCII source'
} finally {
    [Console]::InputEncoding = $previousEncoding
    Remove-Item -LiteralPath $directory -Recurse -Force
}

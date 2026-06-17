# Repair a Windows wheel by bundling extension DLL dependencies.
param(
    [Parameter(Mandatory = $true, Position = 0)]
    [string]$Wheel
)

$ErrorActionPreference = "Stop"

$DestDir = Join-Path $env:GITHUB_WORKSPACE "wheelhouse"
New-Item -ItemType Directory -Force -Path $DestDir | Out-Null

delvewheel repair `
    --add-path (Join-Path $env:VCPKG_ROOT "installed" $env:VCPKG_DEFAULT_TRIPLET "bin") `
    --exclude "python*.dll;vcruntime*.dll;api-ms-win*.dll;ucrtbase.dll" `
    -w $DestDir `
    $Wheel

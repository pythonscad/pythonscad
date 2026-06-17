# Repair a Windows wheel by bundling extension DLL dependencies.
param(
    [Parameter(Mandatory = $true)]
    [string]$Wheel,

    [Parameter(Mandatory = $true)]
    [string]$DestDir
)

$ErrorActionPreference = "Stop"

New-Item -ItemType Directory -Force -Path $DestDir | Out-Null

$VcpkgBin = Join-Path $env:VCPKG_ROOT "installed" $env:VCPKG_DEFAULT_TRIPLET "bin"

delvewheel repair `
    --add-path $VcpkgBin `
    --exclude "python*.dll;vcruntime*.dll;api-ms-win*.dll;ucrtbase.dll" `
    -w $DestDir `
    $Wheel

# Install native build dependencies for Windows wheel builds (MSVC ABI).
# Runs as cibuildwheel [tool.cibuildwheel.windows] before-all.
$ErrorActionPreference = "Stop"

Write-Host "=== install-deps-windows.ps1: bootstrapping vcpkg ==="

# Pin vcpkg to a release tag for reproducible builds (see vcpkg.json builtin-baseline).
$VcpkgVersion = "2025.04.09"
$VcpkgRoot = Join-Path $env:RUNNER_TEMP "vcpkg"

if (-not (Test-Path $VcpkgRoot)) {
    git clone --depth 1 --branch $VcpkgVersion `
        https://github.com/microsoft/vcpkg.git $VcpkgRoot
}

$VcpkgExe = Join-Path $VcpkgRoot "vcpkg.exe"
if (-not (Test-Path $VcpkgExe)) {
    & (Join-Path $VcpkgRoot "bootstrap-vcpkg.bat") -disableMetrics
}

$ManifestDir = Join-Path $env:GITHUB_WORKSPACE "scripts/cibuildwheel"
$Triplet = "x64-windows"

Push-Location $ManifestDir
try {
    & $VcpkgExe install `
        --triplet $Triplet `
        --x-manifest-root $ManifestDir `
        --x-install-root (Join-Path $VcpkgRoot "installed")
} finally {
    Pop-Location
}

$Installed = Join-Path $VcpkgRoot "installed" $Triplet
$PkgConfigDir = Join-Path $Installed "lib" "pkgconfig"
$PkgConfExe = Join-Path $Installed "tools" "pkgconf" "pkgconf.exe"

# Expose vcpkg to setup.py and lex/yacc.
"PYTHONSCAD_VCPKG_ROOT=$VcpkgRoot" | Out-File -FilePath $env:GITHUB_ENV -Encoding utf8 -Append
"VCPKG_ROOT=$VcpkgRoot" | Out-File -FilePath $env:GITHUB_ENV -Encoding utf8 -Append
"VCPKG_DEFAULT_TRIPLET=$Triplet" | Out-File -FilePath $env:GITHUB_ENV -Encoding utf8 -Append
"PKG_CONFIG_PATH=$PkgConfigDir" | Out-File -FilePath $env:GITHUB_ENV -Encoding utf8 -Append
"PKG_CONFIG=$PkgConfExe" | Out-File -FilePath $env:GITHUB_ENV -Encoding utf8 -Append

# winflexbison installs win_bison.exe / win_flex.exe under the triplet tools tree.
$ToolsDir = Join-Path $Installed "tools"
if (Test-Path $ToolsDir) {
    "PATH=$ToolsDir;$env:PATH" | Out-File -FilePath $env:GITHUB_ENV -Encoding utf8 -Append
}

Write-Host "=== install-deps-windows.ps1: done (vcpkg root: $VcpkgRoot) ==="

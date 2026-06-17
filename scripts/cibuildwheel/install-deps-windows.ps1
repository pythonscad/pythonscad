# Install native build dependencies for Windows wheel builds (MSVC ABI).
# Runs as cibuildwheel [tool.cibuildwheel.windows] before-all.
$ErrorActionPreference = "Stop"

Write-Host "=== install-deps-windows.ps1: bootstrapping vcpkg ==="

$ProjectRoot = if ($env:CIBW_PROJECT_DIR) {
    $env:CIBW_PROJECT_DIR
} elseif ($env:GITHUB_WORKSPACE) {
    $env:GITHUB_WORKSPACE
} else {
    (Get-Location).Path
}

# Pin vcpkg to a release tag for reproducible builds (see vcpkg.json builtin-baseline).
$VcpkgVersion = "2025.04.09"
$VcpkgParent = if ($env:RUNNER_TEMP) { $env:RUNNER_TEMP } else { $env:LOCALAPPDATA }
$VcpkgRoot = Join-Path $VcpkgParent "pythonscad-vcpkg"

if (-not (Test-Path $VcpkgRoot)) {
    git clone --depth 1 --branch $VcpkgVersion `
        https://github.com/microsoft/vcpkg.git $VcpkgRoot
}

$VcpkgExe = Join-Path $VcpkgRoot "vcpkg.exe"
if (-not (Test-Path $VcpkgExe)) {
    & (Join-Path $VcpkgRoot "bootstrap-vcpkg.bat") -disableMetrics
}

$ManifestDir = Join-Path $ProjectRoot "scripts/cibuildwheel"
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

# Set env vars for the remainder of this cibuildwheel build (not GITHUB_ENV).
$env:PYTHONSCAD_VCPKG_ROOT = $VcpkgRoot
$env:VCPKG_ROOT = $VcpkgRoot
$env:VCPKG_DEFAULT_TRIPLET = $Triplet
$env:PKG_CONFIG_PATH = $PkgConfigDir
$env:PKG_CONFIG = $PkgConfExe

# winflexbison installs win_bison.exe / win_flex.exe under tools/winflexbison/.
$WinFlexDir = Join-Path $Installed "tools" "winflexbison"
if (Test-Path $WinFlexDir) {
    $env:PATH = "$WinFlexDir;$env:PATH"
}

Write-Host "=== install-deps-windows.ps1: done (vcpkg root: $VcpkgRoot) ==="

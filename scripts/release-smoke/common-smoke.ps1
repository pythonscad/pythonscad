$ErrorActionPreference = 'Stop'

function Write-SmokeLog {
    param([Parameter(Mandatory = $true)][string]$Message)
    Write-Information "==> $Message" -InformationAction Continue
}

function Write-SmokeWarning {
    param([Parameter(Mandatory = $true)][string]$Message)
    Write-Warning $Message
}

function Assert-Command {
    param([Parameter(Mandatory = $true)][string]$Name)
    if (-not (Get-Command $Name -ErrorAction SilentlyContinue)) {
        throw "Required command not found: $Name"
    }
}

function Assert-NonEmptyFile {
    param([Parameter(Mandatory = $true)][string]$Path)
    if (-not (Test-Path -LiteralPath $Path -PathType Leaf)) {
        throw "Expected output file does not exist: $Path"
    }
    if ((Get-Item -LiteralPath $Path).Length -le 0) {
        throw "Expected non-empty output file: $Path"
    }
}

function Write-Utf8Content {
    param(
        [Parameter(Mandatory = $true)][string]$Path,
        [Parameter(Mandatory = $true)][string]$Value
    )
    $encoding = New-Object System.Text.UTF8Encoding -ArgumentList $false
    [System.IO.File]::WriteAllText($Path, $Value, $encoding)
}

function ConvertTo-WindowsArgument {
    param([AllowNull()][string]$Argument)

    if ($null -eq $Argument) {
        return '""'
    }
    if ($Argument -notmatch '[\s"]') {
        return $Argument
    }

    $quoted = '"'
    $backslashes = 0
    foreach ($char in $Argument.ToCharArray()) {
        if ($char -eq '\') {
            $backslashes += 1
        } elseif ($char -eq '"') {
            $quoted += ('\' * (($backslashes * 2) + 1)) + '"'
            $backslashes = 0
        } else {
            if ($backslashes -gt 0) {
                $quoted += '\' * $backslashes
                $backslashes = 0
            }
            $quoted += $char
        }
    }
    if ($backslashes -gt 0) {
        $quoted += '\' * ($backslashes * 2)
    }
    $quoted += '"'
    $quoted
}

function Join-WindowsArgument {
    param([string[]]$Arguments)
    ($Arguments | ForEach-Object { ConvertTo-WindowsArgument -Argument $_ }) -join ' '
}

function Invoke-PythonSCAD {
    param(
        [Parameter(Mandatory = $true)][string]$ExecutablePath,
        [Parameter(Mandatory = $true)][string[]]$Arguments,
        [string]$WorkingDirectory,
        [string]$InputFile,
        [string]$LogFile
    )

    $psi = [System.Diagnostics.ProcessStartInfo]::new()
    $psi.FileName = $ExecutablePath
    $psi.Arguments = Join-WindowsArgument -Arguments $Arguments
    if ($WorkingDirectory) {
        $psi.WorkingDirectory = $WorkingDirectory
    }
    $psi.UseShellExecute = $false
    $psi.Environment['QT_QPA_PLATFORM'] = 'offscreen'
    $psi.RedirectStandardOutput = $true
    $psi.RedirectStandardError = $true
    if ($InputFile) {
        $psi.RedirectStandardInput = $true
    }

    $process = [System.Diagnostics.Process]::Start($psi)
    $stdoutTask = $process.StandardOutput.ReadToEndAsync()
    $stderrTask = $process.StandardError.ReadToEndAsync()
    if ($InputFile) {
        $inputText = Get-Content -LiteralPath $InputFile -Raw
        $process.StandardInput.Write($inputText)
        $process.StandardInput.Close()
    }
    $process.WaitForExit()
    $stdout = $stdoutTask.Result
    $stderr = $stderrTask.Result

    if ($LogFile) {
        $parent = Split-Path -Parent $LogFile
        if ($parent) {
            New-Item -ItemType Directory -Force -Path $parent | Out-Null
        }
        Write-Utf8Content -Path $LogFile -Value ($stdout + $stderr)
    }

    if ($process.ExitCode -ne 0) {
        Write-Error "Smoke command failed with exit code $($process.ExitCode): $ExecutablePath $($Arguments -join ' ')" -ErrorAction Continue
        if ($WorkingDirectory) {
            Write-Error "Working directory: $WorkingDirectory" -ErrorAction Continue
        }
        if ($InputFile) {
            Write-Error "Stdin: $InputFile" -ErrorAction Continue
        }
        if ($LogFile) {
            Write-Error "Captured output: $LogFile" -ErrorAction Continue
            if (Test-Path -LiteralPath $LogFile -PathType Leaf) {
                Write-Error "--- begin command output: $LogFile ---" -ErrorAction Continue
                Get-Content -LiteralPath $LogFile -Raw -ErrorAction SilentlyContinue | Write-Error -ErrorAction Continue
                Write-Error "--- end command output: $LogFile ---" -ErrorAction Continue
            } else {
                Write-Error 'Command produced no captured output.' -ErrorAction Continue
            }
        }
        throw "Smoke command failed with exit code $($process.ExitCode)"
    }
}

function Initialize-SmokeInput {
    param([Parameter(Mandatory = $true)][string]$TestDirectory)

    New-Item -ItemType Directory -Force -Path $TestDirectory | Out-Null
    Write-Utf8Content -Path (Join-Path $TestDirectory 'cube.scad') -Value 'cube(10);'
    Write-Utf8Content -Path (Join-Path $TestDirectory 'cube.py') -Value @'
from pythonscad import *

cube(10).show()
'@
    Write-Utf8Content -Path (Join-Path $TestDirectory 'repl-smoke.py') -Value @'
from pythonscad import *

export(cube(10), "repl-cube.stl")
raise SystemExit(0)
'@
    Write-Utf8Content -Path (Join-Path $TestDirectory 'ipython-smoke.py') -Value @'
from pythonscad import *

export(cube(10), "ipython-cube.stl")
'@
    Write-Utf8Content -Path (Join-Path $TestDirectory 'pyqt6-smoke.py') -Value @'
from PyQt6 import QtCore, QtWidgets, sip
from pythonscad import qapp_ptr

# --repl never starts the GUI, so qapp_ptr() is 0. wrapinstance(0) raises
# and repl() still exits 0, which used to look like a silent packaging miss.
# Importing the bundled bindings is the packaging check; wrap the C++ app
# only when one actually exists.
ptr = qapp_ptr()
if ptr:
    app = sip.wrapinstance(ptr, QtWidgets.QApplication)
    assert app is not None

marker = "PYQT6_PACKAGED_OK " + QtCore.QT_VERSION_STR
print(marker, flush=True)
with open("pyqt6-ok.txt", "w", encoding="utf-8") as fh:
    fh.write(marker + "\n")
raise SystemExit(0)
'@
}

function Invoke-SmokeTest {
    param(
        [Parameter(Mandatory = $true)][string]$ExecutablePath,
        [Parameter(Mandatory = $true)][string]$Label,
        [Parameter(Mandatory = $true)][string]$Workdir,
        [switch]$ExpectPyQt6
    )

    if (-not (Test-Path -LiteralPath $ExecutablePath -PathType Leaf)) {
        throw "Executable not found: $ExecutablePath"
    }

    $safeLabel = $Label -replace '[^A-Za-z0-9_.-]', '_'
    $testdir = Join-Path $Workdir "smoke-$safeLabel"
    Initialize-SmokeInput -TestDirectory $testdir

    Write-SmokeLog "Smoke testing $Label`: startup"
    Invoke-PythonSCAD -ExecutablePath $ExecutablePath `
        -Arguments @('--info') `
        -LogFile (Join-Path $testdir 'info.log')

    Write-SmokeLog "Smoke testing $Label`: OpenSCAD .scad export"
    Invoke-PythonSCAD -ExecutablePath $ExecutablePath `
        -Arguments @('-o', (Join-Path $testdir 'cube-scad.stl'), (Join-Path $testdir 'cube.scad')) `
        -LogFile (Join-Path $testdir 'scad-export.log')
    Assert-NonEmptyFile -Path (Join-Path $testdir 'cube-scad.stl')

    Write-SmokeLog "Smoke testing $Label`: Python CLI export"
    Invoke-PythonSCAD -ExecutablePath $ExecutablePath `
        -Arguments @('--trust-python', '-o', (Join-Path $testdir 'cube-python.stl'), (Join-Path $testdir 'cube.py')) `
        -LogFile (Join-Path $testdir 'python-export.log')
    Assert-NonEmptyFile -Path (Join-Path $testdir 'cube-python.stl')

    Write-SmokeLog "Smoke testing $Label`: basic Python REPL"
    Invoke-PythonSCAD -ExecutablePath $ExecutablePath `
        -Arguments @('--repl') `
        -WorkingDirectory $testdir `
        -InputFile (Join-Path $testdir 'repl-smoke.py') `
        -LogFile (Join-Path $testdir 'repl.log')
    Assert-NonEmptyFile -Path (Join-Path $testdir 'repl-cube.stl')

    Write-SmokeLog "Smoke testing $Label`: IPython"
    Invoke-PythonSCAD -ExecutablePath $ExecutablePath `
        -Arguments @('--ipython', 'ipython-smoke.py') `
        -WorkingDirectory $testdir `
        -LogFile (Join-Path $testdir 'ipython.log')
    Assert-NonEmptyFile -Path (Join-Path $testdir 'ipython-cube.stl')

    if ($ExpectPyQt6) {
        Write-SmokeLog "Smoke testing $Label`: packaged PyQt6"
        $pyqtLog = Join-Path $testdir 'pyqt6.log'
        Invoke-PythonSCAD -ExecutablePath $ExecutablePath `
            -Arguments @('--repl') `
            -WorkingDirectory $testdir `
            -InputFile (Join-Path $testdir 'pyqt6-smoke.py') `
            -LogFile $pyqtLog
        $pyqtOk = Join-Path $testdir 'pyqt6-ok.txt'
        if (-not (Test-Path -LiteralPath $pyqtOk -PathType Leaf) -or ((Get-Item -LiteralPath $pyqtOk).Length -le 0)) {
            Write-Error "--- begin command output: $pyqtLog ---" -ErrorAction Continue
            Get-Content -LiteralPath $pyqtLog -Raw -ErrorAction SilentlyContinue | Write-Error -ErrorAction Continue
            Write-Error "--- end command output: $pyqtLog ---" -ErrorAction Continue
            throw "$Label`: packaged PyQt6 did not emit its success marker"
        }
    }

    Write-SmokeLog "Smoke testing $Label`: OK"
}

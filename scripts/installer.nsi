; Legacy manual NSIS installer script (not used by the CPack-based build).
; The CPack build uses cmake/nsis/NSIS.template.in and CMakeLists.txt instead.
; This script is kept as a reference for direct NSIS builds outside of CMake.

InstallDir ""
; Add cmake/nsis to the include path so bare filenames resolve correctly.
; Run makensis from the repo root or pass -NOCD and set include paths manually.
!addincludedir "..\cmake\nsis"
; Add the bundled UAC plugin directory so UAC.dll is found without a global install.
!addplugindir "..\cmake\nsis\Plugins\x86-unicode"
!include "LogicLib.nsh"
!include "MUI2.nsh"
!include "nsDialogs.nsh"
!include "WinMessages.nsh"
!include "UAC.nsh"
!include "mingw-file-association.nsh"
!include "x64.nsh"
!include "multiuser.nsh"

Name "PythonSCAD"
OutFile "pythonscad_setup.exe"
!define MUI_CUSTOMFUNCTION_GUIINIT GuiInit
RequestExecutionLevel user

; Declare variables defined in multiuser.nsh (AllUsersRadio) and the
; InstallMode variable used across the installer.
; Note: MultiUser.AllUsersRadio is declared inside multiuser.nsh itself.
Var MultiUser.InstallMode

; installer_arch.nsi must define Function .onInit and set architecture-specific
; registry view. See installer32.nsi / installer64.nsi.
!include "installer_arch.nsi"

!insertmacro MUI_PAGE_WELCOME
Page custom MultiUser.InstallModePage_Create MultiUser.InstallModePage_Leave
!insertmacro MUI_PAGE_DIRECTORY
!insertmacro MUI_PAGE_INSTFILES
!insertmacro MUI_PAGE_FINISH
!insertmacro MUI_UNPAGE_CONFIRM
!insertmacro MUI_UNPAGE_INSTFILES
!insertmacro MUI_LANGUAGE "English"

Function GuiInit
  ; Required by UAC plugin for page-based elevation
  !insertmacro UAC_PageElevation_OnGuiInit
FunctionEnd

!define ARP_INS "Software\Microsoft\Windows\CurrentVersion\Uninstall\PythonSCAD"

DirText "Choose a directory to install PythonSCAD"

Section "install"
SetOutPath $INSTDIR
File pythonscad.exe
File /r /x mingw-cross-env examples
File /r /x mingw-cross-env libraries
File /r /x mingw-cross-env fonts
File /r /x mingw-cross-env locale
File /r /x mingw-cross-env color-schemes
File /r /x mingw-cross-env shaders
File /r /x mingw-cross-env templates
${RegisterExtension} "$INSTDIR\pythonscad.exe" ".scad" "PythonSCAD_File"
${RegisterExtension} "$INSTDIR\pythonscad.exe" ".py" "PythonSCAD_Python_File"
CreateDirectory "$SMPROGRAMS\PythonSCAD"
CreateShortCut "$SMPROGRAMS\PythonSCAD\PythonSCAD.lnk" "$INSTDIR\pythonscad.exe"
CreateShortCut "$SMPROGRAMS\PythonSCAD.lnk" "$INSTDIR\pythonscad.exe"
WriteUninstaller "$INSTDIR\Uninstall.exe"
WriteRegStr SHCTX "${ARP_INS}" "DisplayName" "PythonSCAD (remove only)"
WriteRegStr SHCTX "${ARP_INS}" "DisplayVersion" "${VERSION}"
WriteRegStr SHCTX "${ARP_INS}" "Publisher" "The PythonSCAD Developers"
WriteRegStr SHCTX "${ARP_INS}" "URLInfoAbout" "https://pythonscad.org/"
WriteRegStr SHCTX "${ARP_INS}" "URLUpdateInfo" "https://pythonscad.org/downloads.html"
WriteRegStr SHCTX "${ARP_INS}" "HelpLink" "https://pythonscad.org/"
WriteRegStr SHCTX "${ARP_INS}" "UninstallString" "$\"$INSTDIR\Uninstall.exe$\""
; Store install mode so the uninstaller can determine scope and elevate if needed
${If} $MultiUser.InstallMode = 1
  WriteRegStr HKLM "${ARP_INS}" "InstallMode" "AllUsers"
${Else}
  WriteRegStr HKCU "${ARP_INS}" "InstallMode" "CurrentUser"
${EndIf}
SectionEnd

Section "Uninstall"
Call un.MultiUser.Init
${UnRegisterExtension} ".scad" "PythonSCAD_File"
${UnRegisterExtension} ".py" "PythonSCAD_Python_File"
Delete "$INSTDIR\Uninstall.exe"
Delete "$INSTDIR\pythonscad.exe"
Delete "$SMPROGRAMS\PythonSCAD\PythonSCAD.lnk"
RMDir "$SMPROGRAMS\PythonSCAD"
Delete "$SMPROGRAMS\PythonSCAD.lnk"
DeleteRegKey SHCTX "SOFTWARE\Microsoft\Windows\CurrentVersion\Uninstall\PythonSCAD"
RMDir /r "$INSTDIR\fonts"
RMDir /r "$INSTDIR\color-schemes"
RMDir /r "$INSTDIR\templates"
RMDir /r "$INSTDIR\examples"
RMDir /r "$INSTDIR\libraries\mcad"
RMDir /r "$INSTDIR\locale"
RMDir /r "$INSTDIR\shaders"
Delete "$INSTDIR\libraries\boxes.scad"
Delete "$INSTDIR\libraries\shapes.scad"
RMDir "$INSTDIR\libraries"
RMDir "$INSTDIR"
SectionEnd

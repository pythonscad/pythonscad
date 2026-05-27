; Legacy manual NSIS installer script (not used by the CPack-based build).
; The CPack build uses cmake/nsis/NSIS.template.in and CMakeLists.txt instead.
; This script is kept as a reference for direct NSIS builds outside of CMake.
;
; To build for 64-bit, run:
;   makensis /DARCH=x64 /I..\cmake\nsis scripts\installer.nsi
; To build for 32-bit, run:
;   makensis /DARCH=x86 /I..\cmake\nsis scripts\installer.nsi

InstallDir ""
; Add cmake/nsis to the include path so bare filenames resolve correctly.
!addincludedir "..\cmake\nsis"
; Add the bundled UAC plugin directory so UAC.dll is found without a global install.
!addplugindir "..\cmake\nsis\Plugins\x86-unicode"

; MUI_CUSTOMFUNCTION_GUIINIT must be defined before !include MUI2.nsh
!define MUI_CUSTOMFUNCTION_GUIINIT GuiInit

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
RequestExecutionLevel user

; $MultiUser.InstallMode declared here; MultiUser.AllUsersRadio is in multiuser.nsh
Var MultiUser.InstallMode

; Include architecture-specific .onInit (provides SetRegView and UAC init).
; Pass /DARCH=x64 or /DARCH=x86 to makensis to select the right file.
!ifdef ARCH
  !if "${ARCH}" == "x64"
    !include "installer64.nsi"
  !else
    !include "installer32.nsi"
  !endif
!else
  !warning "No /DARCH flag specified; defaulting to 32-bit init. Pass /DARCH=x64 for 64-bit."
  !include "installer32.nsi"
!endif

!insertmacro MUI_PAGE_WELCOME
Page custom MultiUser.InstallModePage_Create MultiUser.InstallModePage_Leave
!insertmacro MUI_PAGE_DIRECTORY
!insertmacro MUI_PAGE_INSTFILES
!insertmacro MUI_PAGE_FINISH
!insertmacro MUI_UNPAGE_CONFIRM
!insertmacro MUI_UNPAGE_INSTFILES
!insertmacro MUI_LANGUAGE "English"

Function GuiInit
  !insertmacro UAC_PageElevation_OnGuiInit
FunctionEnd

!define ARP_INS "Software\Microsoft\Windows\CurrentVersion\Uninstall\PythonSCAD"

DirText "Choose a directory to install PythonSCAD"

Section "install"
SetOutPath $INSTDIR
File pythonscad.exe
File pythonscad.com
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

Function un.onInit
  ; Match the 64-bit registry view used at install time so uninstall keys are found.
  ${If} ${RunningX64}
    SetRegView 64
  ${EndIf}
  Call un.MultiUser.Init
  ; Elevate for all-users uninstall if not already admin
  ${If} $MultiUser.InstallMode = 1
  ${AndIfNot} ${UAC_IsAdmin}
    !insertmacro UAC_RunElevated
    ${Switch} $0
    ${Case} 0
      ${IfThen} $1 = 1 ${|} Quit ${|}   ; outer instance, inner did the work
      ${IfThen} $3 <> 0 ${|} ${Break} ${|}  ; we are admin, continue
      MessageBox MB_ICONSTOP|MB_TOPMOST "Administrator privileges are required to uninstall for all users."
      Quit
    ${Case} 1223
      Quit   ; User cancelled UAC
    ${Case} 1062
      MessageBox MB_ICONSTOP|MB_TOPMOST "Logon service not running, cannot elevate."
      Quit
    ${Default}
      MessageBox MB_ICONSTOP|MB_TOPMOST "Unable to elevate, error $0"
      Quit
    ${EndSwitch}
  ${EndIf}
FunctionEnd

Section "Uninstall"
${UnRegisterExtension} ".scad" "PythonSCAD_File"
${UnRegisterExtension} ".py" "PythonSCAD_Python_File"
Delete "$INSTDIR\Uninstall.exe"
Delete "$INSTDIR\pythonscad.exe"
Delete "$INSTDIR\pythonscad.com"
Delete "$SMPROGRAMS\PythonSCAD\PythonSCAD.lnk"
RMDir "$SMPROGRAMS\PythonSCAD"
Delete "$SMPROGRAMS\PythonSCAD.lnk"
DeleteRegKey SHCTX "SOFTWARE\Microsoft\Windows\CurrentVersion\Uninstall\PythonSCAD"
RMDir /r "$INSTDIR\fonts"
RMDir /r "$INSTDIR\color-schemes"
RMDir /r "$INSTDIR\templates"
RMDir /r "$INSTDIR\examples"
RMDir /r "$INSTDIR\libraries"
RMDir /r "$INSTDIR\locale"
RMDir /r "$INSTDIR\shaders"
RMDir "$INSTDIR"
SectionEnd

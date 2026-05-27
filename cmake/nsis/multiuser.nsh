; Multiuser installation support for PythonSCAD
; Based on the NSIS UAC plugin DualMode example (zlib license)
;
; $MultiUser.InstallMode: 0 = current user only, 1 = all users
; (Var MultiUser.InstallMode is declared in CMakeLists.txt CPACK_NSIS_DEFINES
;  so it appears at the top level before this file is included.)

!include "WinMessages.nsh"

!ifndef BCM_SETSHIELD
!define BCM_SETSHIELD 0x0000160C
!endif

; Default install directories
!ifndef MULTIUSER_INSTALLDIR_USER
  !define MULTIUSER_INSTALLDIR_USER "$LocalAppData\PythonSCAD"
!endif
; Default all-users install dir. Must be defined BEFORE !include "multiuser.nsh"
; to override (e.g. installer64.nsi defines it to use $PROGRAMFILES64).
; $PROGRAMFILES resolves to Program Files (x86) on x64 Windows when the NSIS
; process is 32-bit, so x64 builds must define MULTIUSER_INSTALLDIR_ALLUSERS
; to $PROGRAMFILES64\PythonSCAD before including this file.
!ifndef MULTIUSER_INSTALLDIR_ALLUSERS
  !define MULTIUSER_INSTALLDIR_ALLUSERS "$PROGRAMFILES\PythonSCAD"
!endif

; Dedicated variable for the "all users" radio button HWND, avoiding
; cross-function stack passing which is fragile and error-prone.
Var MultiUser.AllUsersRadio

; --- Apply shell context based on current mode ---
; Called at the start of the install Section (via CPACK_NSIS_EXTRA_PREINSTALL_COMMANDS)
; and from the install-type page Leave function.
Function MultiUser.SetContext
  ${If} $MultiUser.InstallMode = 0
    SetShellVarContext current
  ${Else}
    SetShellVarContext all
  ${EndIf}
FunctionEnd

; --- Radio button click callback: update UAC shield on Next button ---
; Show shield only when "For all users" is selected AND we are not already admin.
Function MultiUser.InstallModePage_OnClick
  Pop $1
  nsDialogs::GetUserData $1
  Pop $1
  ; If already elevated, never show the shield — no further UAC prompt will occur.
  ${If} ${UAC_IsAdmin}
    StrCpy $1 0
  ${EndIf}
  GetDlgItem $0 $HWNDParent 1
  SendMessage $0 ${BCM_SETSHIELD} 0 $1
FunctionEnd

; --- Remove shield when user goes Back ---
Function MultiUser.RemoveShield
  GetDlgItem $0 $HWNDParent 1
  SendMessage $0 ${BCM_SETSHIELD} 0 0
FunctionEnd

; --- Install mode selection page ---
Function MultiUser.InstallModePage_Create
  !insertmacro MUI_HEADER_TEXT "Installation Type" "Choose how PythonSCAD should be installed."
  GetFunctionAddress $8 MultiUser.InstallModePage_OnClick
  nsDialogs::Create /NOUNLOAD 1018
  Pop $9
  ${NSD_OnBack} MultiUser.RemoveShield

  ${NSD_CreateLabel} 0 10u 100% 20u "Select whether to install PythonSCAD for yourself only or for all users on this computer."
  Pop $0

  ; "Just for me" radio button
  ${NSD_CreateRadioButton} 10u 38u 90% 15u "Just for me (no administrator rights required)"
  Pop $0
  nsDialogs::OnClick $0 $8
  nsDialogs::SetUserData $0 0
  ${IfThen} $MultiUser.InstallMode = 0 ${|} SendMessage $0 ${BM_CLICK} 0 0 ${|}

  ; "For all users" radio button — store HWND in dedicated variable
  ${NSD_CreateRadioButton} 10u 58u 90% 15u "For all users on this computer (requires administrator rights)"
  Pop $MultiUser.AllUsersRadio
  nsDialogs::OnClick $MultiUser.AllUsersRadio $8
  nsDialogs::SetUserData $MultiUser.AllUsersRadio 1
  ${IfThen} $MultiUser.InstallMode = 1 ${|} SendMessage $MultiUser.AllUsersRadio ${BM_CLICK} 0 0 ${|}

  nsDialogs::Show
FunctionEnd

; --- Leave callback: set mode, update InstDir, and elevate if needed ---
Function MultiUser.InstallModePage_Leave
  ${NSD_GetState} $MultiUser.AllUsersRadio $9

  ${If} $9 = 0
    ; Just for me
    StrCpy $MultiUser.InstallMode 0
    StrCpy $INSTDIR "${MULTIUSER_INSTALLDIR_USER}"
    Call MultiUser.SetContext
  ${Else}
    ; All users — elevate if not already admin
    StrCpy $MultiUser.InstallMode 1
    StrCpy $INSTDIR "${MULTIUSER_INSTALLDIR_ALLUSERS}"
    Call MultiUser.SetContext
    ${If} ${UAC_IsAdmin}
      ; Already elevated: clear the shield from the Next button since no UAC is needed.
      Call MultiUser.RemoveShield
    ${Else}
      GetDlgItem $9 $HWNDParent 1
      System::Call 'user32::GetFocus()i.r8'
      EnableWindow $9 0
      !insertmacro UAC_PageElevation_RunElevated
      EnableWindow $9 1
      System::Call 'user32::SetFocus(ir8)'
      ${If} $2 = 0x666666
        MessageBox MB_ICONEXCLAMATION "Administrator privileges are required to install for all users.$\nPlease select 'Just for me' or log in with an administrator account."
        Abort
      ${ElseIf} $0 = 1223
        Abort   ; User cancelled UAC dialog
      ${Else}
        ${If} $0 <> 0
          ${If} $0 = 1062
            MessageBox MB_ICONSTOP "Unable to elevate: Secondary Logon service is not running."
          ${Else}
            MessageBox MB_ICONSTOP "Unable to elevate, error code: $0"
          ${EndIf}
          Abort
        ${EndIf}
      ${EndIf}
      Quit  ; Outer instance exits; elevated inner instance continues
    ${EndIf}
  ${EndIf}
FunctionEnd

; --- Uninstaller: read stored install mode and set shell context ---
; Scope is determined by matching the running uninstaller path ($EXEPATH) against
; the UninstallString stored in HKLM (all-users) and HKCU (per-user). This
; correctly disambiguates scope when both install types exist on the same machine.
; Note: HKCU is the current user's hive only — per-user installs for other accounts
; are not visible and cannot be uninstalled by this uninstaller.
!define ARP_UN "Software\Microsoft\Windows\CurrentVersion\Uninstall\PythonSCAD"

Function un.MultiUser.Init
  ; Determine install scope by matching the running uninstaller path against the
  ; UninstallString stored in each hive. This correctly handles the case where
  ; both a per-user and an all-users install exist simultaneously, avoiding the
  ; HKLM-first ambiguity that would misclassify a per-user uninstall.
  ReadRegStr $0 HKLM "${ARP_UN}" "UninstallString"
  ${If} $0 != ""
    ; Strip surrounding quotes if present (leading then trailing, independently).
    StrCpy $1 $0 1 0
    ${If} $1 == "$\""
      StrCpy $0 $0 "" 1   ; remove leading quote
    ${EndIf}
    StrCpy $1 $0 1 -1
    ${If} $1 == "$\""
      StrCpy $0 $0 -1     ; remove trailing quote
    ${EndIf}
    ${If} $0 == "$EXEPATH"
      StrCpy $MultiUser.InstallMode 1
      SetShellVarContext all
      Return
    ${EndIf}
  ${EndIf}
  ; Check HKCU for per-user install
  ReadRegStr $0 HKCU "${ARP_UN}" "UninstallString"
  ${If} $0 != ""
    ; Strip surrounding quotes if present (leading then trailing, independently).
    StrCpy $1 $0 1 0
    ${If} $1 == "$\""
      StrCpy $0 $0 "" 1   ; remove leading quote
    ${EndIf}
    StrCpy $1 $0 1 -1
    ${If} $1 == "$\""
      StrCpy $0 $0 -1     ; remove trailing quote
    ${EndIf}
    ${If} $0 == "$EXEPATH"
      StrCpy $MultiUser.InstallMode 0
      SetShellVarContext current
      Return
    ${EndIf}
  ${EndIf}
  ; Default to current user if no matching key found
  StrCpy $MultiUser.InstallMode 0
  SetShellVarContext current
FunctionEnd

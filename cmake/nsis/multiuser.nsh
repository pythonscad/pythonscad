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
; Use $PROGRAMFILES (not $PROGRAMFILES64) so Windows selects the correct
; Program Files directory for the installer's bitness automatically.
; Override this from CMakeLists.txt via CPACK_NSIS_DEFINES if needed.
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
Function MultiUser.InstallModePage_OnClick
  Pop $1
  nsDialogs::GetUserData $1
  Pop $1
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
    ${IfNot} ${UAC_IsAdmin}
      GetDlgItem $9 $HWNDParent 1
      System::Call 'user32::GetFocus()i.s'
      EnableWindow $9 0
      !insertmacro UAC_PageElevation_RunElevated
      EnableWindow $9 1
      System::Call 'user32::SetFocus(is)'
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
; InstallMode is stored under the ARP uninstall key in the appropriate hive
; (HKLM for all-users, HKCU for per-user). Reading from both hives in the
; right order ensures correct scope even when called from a different user account.
!define ARP_UN "Software\Microsoft\Windows\CurrentVersion\Uninstall\PythonSCAD"

Function un.MultiUser.Init
  ; Try HKLM first — written only for all-users installs, readable by any account
  ReadRegStr $0 HKLM "${ARP_UN}" "InstallMode"
  ${If} $0 == "AllUsers"
    StrCpy $MultiUser.InstallMode 1
    SetShellVarContext all
    Return
  ${EndIf}
  ; Fall back to HKCU for per-user installs
  ReadRegStr $0 HKCU "${ARP_UN}" "InstallMode"
  ${If} $0 == "CurrentUser"
    StrCpy $MultiUser.InstallMode 0
    SetShellVarContext current
    Return
  ${EndIf}
  ; Default to current user if no key found
  StrCpy $MultiUser.InstallMode 0
  SetShellVarContext current
FunctionEnd

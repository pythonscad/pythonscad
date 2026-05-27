; Multiuser installation support for PythonSCAD
; Based on the NSIS UAC plugin DualMode example (zlib license)
;
; $MultiUser.InstallMode: 0 = current user only, 1 = all users
; (Var MultiUser.InstallMode is declared in CMakeLists.txt CPACK_NSIS_DEFINES
;  so it appears at the top level before this file is included.)

!ifndef BCM_SETSHIELD
!define BCM_SETSHIELD 0x0000160C
!endif

; Default install directories
!ifndef MULTIUSER_INSTALLDIR_USER
  !define MULTIUSER_INSTALLDIR_USER "$LocalAppData\PythonSCAD"
!endif
!ifndef MULTIUSER_INSTALLDIR_ALLUSERS
  !define MULTIUSER_INSTALLDIR_ALLUSERS "$ProgramFiles64\PythonSCAD"
!endif

; --- Apply shell context and default InstDir based on current mode ---
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

  ; "For all users" radio button (hwnd stored in $2 for Leave to read)
  ${NSD_CreateRadioButton} 10u 58u 90% 15u "For all users on this computer (requires administrator rights)"
  Pop $2
  nsDialogs::OnClick $2 $8
  nsDialogs::SetUserData $2 1
  ${IfThen} $MultiUser.InstallMode = 1 ${|} SendMessage $2 ${BM_CLICK} 0 0 ${|}

  Push $2   ; keep "all users" hwnd on stack for Leave
  nsDialogs::Show
  Pop $2
FunctionEnd

; --- Leave callback: set mode, update InstDir, and elevate if needed ---
Function MultiUser.InstallModePage_Leave
  ; Stack top has the "all users" radio hwnd (pushed at end of _Create)
  Pop $0
  Push $0
  ${NSD_GetState} $0 $9

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
      System::Call user32::GetFocus()i.s
      EnableWindow $9 0
      !insertmacro UAC_PageElevation_RunElevated
      EnableWindow $9 1
      System::Call user32::SetFocus(is)
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
; Called from un.onInit via CPACK_NSIS_EXTRA_UNINSTALL_COMMANDS
Function un.MultiUser.Init
  ReadRegStr $0 HKCU "Software\PythonSCAD" "InstallMode"
  ${If} $0 == "AllUsers"
    StrCpy $MultiUser.InstallMode 1
    SetShellVarContext all
  ${Else}
    StrCpy $MultiUser.InstallMode 0
    SetShellVarContext current
  ${EndIf}
FunctionEnd

Function .onInit
  ${If} ${RunningX64}
    SetRegView 64
  ${Else}
    MessageBox MB_OK "This is the 64-bit PythonSCAD installer. Your system appears to be 32-bit.$\nPlease download the 32-bit installer instead."
    Quit
  ${EndIf}
  !insertmacro UAC_PageElevation_OnInit
  StrCpy $MultiUser.InstallMode 0
  ${IfThen} ${UAC_IsAdmin} ${|} StrCpy $MultiUser.InstallMode 1 ${|}
FunctionEnd

Function .onInit
  !insertmacro UAC_PageElevation_OnInit
  StrCpy $MultiUser.InstallMode 0
  ${IfThen} ${UAC_IsAdmin} ${|} StrCpy $MultiUser.InstallMode 1 ${|}
FunctionEnd

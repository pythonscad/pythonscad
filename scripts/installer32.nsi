Function .onInit
  !insertmacro UAC_PageElevation_OnInit
  ${IfThen} ${UAC_IsAdmin} ${|} StrCpy $MultiUser.InstallMode 1 ${|}
FunctionEnd

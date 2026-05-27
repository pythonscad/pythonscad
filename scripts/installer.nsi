InstallDir ""
!include "LogicLib.nsh"
!include "MUI2.nsh"
!include "nsDialogs.nsh"
!include "UAC.nsh"
!include "mingw-file-association.nsh"
!include "x64.nsh"
!include "multiuser.nsh"

Name "PythonSCAD"
OutFile "pythonscad_setup.exe"
!define MUI_CUSTOMFUNCTION_GUIINIT GuiInit
RequestExecutionLevel user

!include "installer_arch.nsi"

page custom MultiUser.InstallModePage_Create MultiUser.InstallModePage_Leave

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
WriteUninstaller $INSTDIR\Uninstall.exe
WriteRegStr SHCTX "Software\Microsoft\Windows\CurrentVersion\Uninstall\PythonSCAD" "DisplayName" "PythonSCAD (remove only)"
WriteRegStr SHCTX "Software\Microsoft\Windows\CurrentVersion\Uninstall\PythonSCAD" "DisplayVersion" "${VERSION}"
WriteRegStr SHCTX "Software\Microsoft\Windows\CurrentVersion\Uninstall\PythonSCAD" "Publisher" "The PythonSCAD Developers"
WriteRegStr SHCTX "Software\Microsoft\Windows\CurrentVersion\Uninstall\PythonSCAD" "URLInfoAbout" "https://pythonscad.org/"
WriteRegStr SHCTX "Software\Microsoft\Windows\CurrentVersion\Uninstall\PythonSCAD" "URLUpdateInfo" "https://pythonscad.org/downloads.html"
WriteRegStr SHCTX "Software\Microsoft\Windows\CurrentVersion\Uninstall\PythonSCAD" "HelpLink" "https://pythonscad.org/"
WriteRegStr SHCTX "Software\Microsoft\Windows\CurrentVersion\Uninstall\PythonSCAD" "UninstallString" "$INSTDIR\Uninstall.exe"
SectionEnd

Section "Uninstall"
Call un.MultiUser.Init
${UnRegisterExtension} ".scad" "PythonSCAD_File"
${UnRegisterExtension} ".py" "PythonSCAD_Python_File"
Delete $INSTDIR\Uninstall.exe
Delete $INSTDIR\pythonscad.exe
Delete "$SMPROGRAMS\PythonSCAD\PythonSCAD.lnk"
RMDir "$SMPROGRAMS\PythonSCAD"
Delete "$SMPROGRAMS\PythonSCAD.lnk"
DeleteRegKey SHCTX "SOFTWARE\Microsoft\Windows\CurrentVersion\Uninstall\PythonSCAD"
RMDir /r $INSTDIR\fonts
RMDir /r $INSTDIR\color-schemes
RMDir /r $INSTDIR\templates
RMDir /r $INSTDIR\examples
RMDir /r $INSTDIR\libraries\mcad
RMDir /r $INSTDIR\locale
RMDir /r $INSTDIR\shaders
Delete $INSTDIR\libraries\boxes.scad
Delete $INSTDIR\libraries\shapes.scad
RMDir $INSTDIR\libraries
RMDir $INSTDIR
SectionEnd

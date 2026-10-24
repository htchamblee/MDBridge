; Inno Setup script of the MDBridge Windows installer.
;
;   ISCC.exe /DAppVersion=0.1.0 /DSourceDir=<install tree> /DOutputDir=<dist> ^
;            /DOutputName=MDBridge-0.1.0-windows-x64-setup /DIconFile=<mdbridge.ico> packaging\windows\mdbridge.iss
;
; packaging\windows\build_windows.ps1 calls it with these values; <install tree> is what
; `cmake --install` wrote (bin\ with the programs, Qt and the MSVC runtime; licences beside it).
; Installs per user by default (no administrator rights needed); "Install for all users" is offered
; in the setup dialog.

#ifndef AppVersion
  #define AppVersion "0.0.0"
#endif
#ifndef SourceDir
  #error SourceDir (the install tree of the build) must be given with /DSourceDir=...
#endif
#ifndef OutputDir
  #define OutputDir "."
#endif
#ifndef OutputName
  #define OutputName "MDBridge-" + AppVersion + "-windows-x64-setup"
#endif

[Setup]
AppId={{5D17C581-B109-4915-BCF0-81DA7AEE87DA}
AppName=MDBridge
AppVersion={#AppVersion}
AppVerName=MDBridge {#AppVersion}
AppPublisher=MDBridge
AppPublisherURL=https://github.com/htchamblee/MDBridge
AppSupportURL=https://github.com/htchamblee/MDBridge/issues
AppComments=Analysis and visualisation of MD trajectories of RNA systems
AppCopyright=Copyright (c) 2026 the MDBridge author. All rights reserved.
DefaultDirName={autopf}\MDBridge
DefaultGroupName=MDBridge
DisableProgramGroupPage=yes
LicenseFile={#SourceDir}\LICENSE.txt
OutputDir={#OutputDir}
OutputBaseFilename={#OutputName}
#ifdef IconFile
SetupIconFile={#IconFile}
#endif
UninstallDisplayIcon={app}\bin\mdbridge-gui.exe
Compression=lzma2/max
SolidCompression=yes
WizardStyle=modern
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
MinVersion=10.0
PrivilegesRequired=lowest
PrivilegesRequiredOverridesAllowed=dialog
ChangesEnvironment=yes

[Tasks]
Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}"
Name: "addtopath"; Description: "Add the mdbridge command line program to PATH"; GroupDescription: "Command line:"; Flags: unchecked

[Files]
Source: "{#SourceDir}\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs

[Icons]
Name: "{group}\MDBridge"; Filename: "{app}\bin\mdbridge-gui.exe"
Name: "{group}\MDBridge command prompt"; Filename: "{cmd}"; Parameters: "/k ""{app}\bin\mdbridge.exe"" --help"; \
    WorkingDir: "{%USERPROFILE}"
Name: "{group}\Licences"; Filename: "{app}\THIRD_PARTY_NOTICES.txt"
Name: "{group}\{cm:UninstallProgram,MDBridge}"; Filename: "{uninstallexe}"
Name: "{autodesktop}\MDBridge"; Filename: "{app}\bin\mdbridge-gui.exe"; Tasks: desktopicon

[Registry]
; per-user or all-users PATH entry, depending on the install mode
Root: HKCU; Subkey: "Environment"; ValueType: expandsz; ValueName: "Path"; ValueData: "{olddata};{app}\bin"; \
    Tasks: addtopath; Check: (not IsAdminInstallMode) and NeedsAddPath(ExpandConstant('{app}\bin'))
Root: HKLM; Subkey: "SYSTEM\CurrentControlSet\Control\Session Manager\Environment"; ValueType: expandsz; \
    ValueName: "Path"; ValueData: "{olddata};{app}\bin"; \
    Tasks: addtopath; Check: IsAdminInstallMode and NeedsAddPath(ExpandConstant('{app}\bin'))

[Run]
Filename: "{app}\bin\mdbridge-gui.exe"; Description: "{cm:LaunchProgram,MDBridge}"; Flags: nowait postinstall skipifsilent

[Code]
function EnvironmentKey(var Root: Integer): string;
begin
  if IsAdminInstallMode then begin
    Root := HKEY_LOCAL_MACHINE;
    Result := 'SYSTEM\CurrentControlSet\Control\Session Manager\Environment';
  end else begin
    Root := HKEY_CURRENT_USER;
    Result := 'Environment';
  end;
end;

function NeedsAddPath(Dir: string): Boolean;
var
  Root: Integer;
  Key, Paths: string;
begin
  Key := EnvironmentKey(Root);
  if not RegQueryStringValue(Root, Key, 'Path', Paths) then begin
    Result := True;
    exit;
  end;
  Result := Pos(';' + Uppercase(Dir) + ';', ';' + Uppercase(Paths) + ';') = 0;
end;

procedure RemovePath(Dir: string);
var
  Root: Integer;
  Key, Paths: string;
  P: Integer;
begin
  Key := EnvironmentKey(Root);
  if not RegQueryStringValue(Root, Key, 'Path', Paths) then
    exit;
  P := Pos(';' + Uppercase(Dir), Uppercase(Paths));
  if P > 0 then begin
    Delete(Paths, P, Length(Dir) + 1);
    RegWriteExpandStringValue(Root, Key, 'Path', Paths);
  end;
end;

procedure CurUninstallStepChanged(CurUninstallStep: TUninstallStep);
begin
  if CurUninstallStep = usPostUninstall then
    RemovePath(ExpandConstant('{app}\bin'));
end;

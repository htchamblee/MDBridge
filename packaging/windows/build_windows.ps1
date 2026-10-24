# Build MDBridge for Windows x64: programs, tests, zip and installer.
#
#   packaging\windows\build_windows.ps1 -QtDir C:\Qt\6.11.3\msvc2022_64 [-Out dist] [-Work build\windows] [-SkipTests]
#
# Needs Visual Studio 2022 with the C++ workload, CMake 3.24+, Qt 6 for MSVC 2022 (-QtDir), and
# Inno Setup 6 for the installer (skipped with a warning when ISCC.exe is not found).
# Result: OUT\MDBridge-<version>-windows-x64.zip and OUT\MDBridge-<version>-windows-x64-setup.exe;
# the unpacked programs are in WORK\install.
param(
    [Parameter(Mandatory = $true)][string]$QtDir,
    [string]$Out = "dist",
    [string]$Work = "build\windows",
    [switch]$SkipTests
)
$ErrorActionPreference = "Stop"
$Root = (Resolve-Path "$PSScriptRoot\..\..").Path
$Build = "$Work\build"
$Install = "$Work\install"

function Invoke-Checked {
    & $args[0] @($args | Select-Object -Skip 1)
    if ($LASTEXITCODE -ne 0) { throw "failed (exit $LASTEXITCODE): $($args -join ' ')" }
}

Write-Output "== configure and build (Qt at $QtDir)"
Invoke-Checked cmake -S $Root -B $Build -G "Visual Studio 17 2022" -A x64 "-DCMAKE_PREFIX_PATH=$QtDir" `
    -DMDBRIDGE_WARNINGS_AS_ERRORS=ON
Invoke-Checked cmake --build $Build --config Release --parallel

if (-not $SkipTests) {
    Write-Output "== tests"
    $env:PATH = "$QtDir\bin;$env:PATH"  # the test programs use Qt's DLLs from the Qt folder
    Invoke-Checked ctest --test-dir $Build -C Release --output-on-failure
}

Write-Output "== install (windeployqt, MSVC runtime, licences)"
if (Test-Path $Install) { Remove-Item -Recurse -Force $Install }
Invoke-Checked cmake --install $Build --config Release --prefix $Install
$version = (Select-String -Path "$Build\CMakeCache.txt" -Pattern '^CMAKE_PROJECT_VERSION:STATIC=(.+)$').Matches[0].Groups[1].Value
$name = "MDBridge-$version-windows-x64"
New-Item -ItemType Directory -Force $Out | Out-Null
$Out = (Resolve-Path $Out).Path

Write-Output "== zip"
$zipRoot = "$Work\zip\$name"
if (Test-Path "$Work\zip") { Remove-Item -Recurse -Force "$Work\zip" }
New-Item -ItemType Directory -Force $zipRoot | Out-Null
Copy-Item -Recurse "$Install\*" $zipRoot
Compress-Archive -Path $zipRoot -DestinationPath "$Out\$name.zip" -Force

Write-Output "== installer"
$iscc = @("${env:ProgramFiles(x86)}\Inno Setup 6\ISCC.exe", "$env:ProgramFiles\Inno Setup 6\ISCC.exe") |
    Where-Object { Test-Path $_ } | Select-Object -First 1
if ($iscc) {
    Invoke-Checked $iscc /Q "/DAppVersion=$version" "/DSourceDir=$((Resolve-Path $Install).Path)" "/DOutputDir=$Out" `
        "/DOutputName=$name-setup" "/DIconFile=$Root\data\icons\mdbridge.ico" "$Root\packaging\windows\mdbridge.iss"
} else {
    Write-Warning "Inno Setup 6 (ISCC.exe) not found: no installer"
}
Get-ChildItem $Out | Format-Table Name, Length

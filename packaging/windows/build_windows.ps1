# Build MDBridge for Windows x64: programs, tests, zip and installer.
#
#   packaging\windows\build_windows.ps1 -QtDir C:\Qt\6.11.3\msvc2022_64 [-Out dist] [-Work build\windows] [-SkipTests]
#                                         [-FFmpegRoot build\ffmpeg\prefix]
#
# Needs Visual Studio 2022 with the C++ workload, CMake 3.24+, Qt 6 for MSVC 2022 (-QtDir), and
# Inno Setup 6 for the installer (skipped with a warning when ISCC.exe is not found). -FFmpegRoot:
# the LGPL FFmpeg libraries of packaging\ffmpeg\windows\build_ffmpeg.ps1 (required then; without
# it, video encoding is used when CMake finds the libraries, else left out).
# Result: OUT\MDBridge-<version>-windows-x64.zip and OUT\MDBridge-<version>-windows-x64-setup.exe;
# the unpacked programs are in WORK\install.
param(
    [Parameter(Mandatory = $true)][string]$QtDir,
    [string]$Out = "dist",
    [string]$Work = "build\windows",
    [switch]$SkipTests,
    [string]$FFmpegRoot = ""
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
$ffmpegArgs = @()
if ($FFmpegRoot) {
    $ffmpegArgs = @("-DMDBRIDGE_FFMPEG=ON", "-DMDBRIDGE_FFMPEG_ROOT=$((Resolve-Path $FFmpegRoot).Path)")
}
Invoke-Checked cmake -S $Root -B $Build -G "Visual Studio 17 2022" -A x64 "-DCMAKE_PREFIX_PATH=$QtDir" `
    -DMDBRIDGE_WARNINGS_AS_ERRORS=ON @ffmpegArgs
Invoke-Checked cmake --build $Build --config Release --parallel

if (-not $SkipTests) {
    Write-Output "== tests"
    $env:PATH = "$QtDir\bin;$env:PATH"  # the test programs use Qt's DLLs from the Qt folder
    Invoke-Checked ctest --test-dir $Build -C Release --output-on-failure
}

Write-Output "== install (windeployqt, MSVC runtime, licences)"
if (Test-Path $Install) { Remove-Item -Recurse -Force $Install }
Invoke-Checked cmake --install $Build --config Release --prefix $Install
# JPEG figures (charts / preview as .jpg) need Qt's JPEG image plugin, which windeployqt deploys with
# the other image formats; the zip and the installer take the install folder as it is
if (-not (Get-ChildItem -Path $Install -Recurse -Filter qjpeg.dll | Where-Object { $_.Directory.Name -eq "imageformats" })) {
    throw "Qt's JPEG image plugin (imageformats\qjpeg.dll) is missing from the install: JPEG charts would fail"
}
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

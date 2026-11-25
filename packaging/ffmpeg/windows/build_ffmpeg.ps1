# Build MDBridge's video encoding libraries on Windows with MSVC: SVT-AV1 and an LGPL-only FFmpeg,
# through vcpkg's ports (which bring FFmpeg's MSYS2 build shell and NASM themselves).
#
#   packaging\ffmpeg\windows\build_ffmpeg.ps1 [-Prefix build\ffmpeg\prefix] [-Work build\ffmpeg\work]
#
# Result, found by CMake with -DMDBRIDGE_FFMPEG_ROOT=<Prefix> (the layout of the Linux script):
#   Prefix\include, Prefix\lib (import libraries), Prefix\bin (avcodec-*.dll, avformat-*.dll,
#   avutil-*.dll, swscale-*.dll, SvtAv1Enc.dll: bundled next to mdbridge.exe)
#   Prefix\share\mdbridge-ffmpeg\   BUILDINFO.txt and the licences of both projects
#   Prefix\..\ffmpeg-<version>-vcpkg-source.zip
#                                   the corresponding source of the FFmpeg DLLs (LGPL 2.1 section 6):
#                                   FFmpeg's source archive as vcpkg downloaded it, vcpkg's ffmpeg port
#                                   (portfile and patches), this script and versions.env
#
# Versions: packaging\ffmpeg\versions.env. vcpkg is cloned at the pinned release tag and commit; its
# ffmpeg and svt-av1 ports must be the pinned FFmpeg and SVT-AV1 versions (checked). vcpkg checks
# the hashes of the source archives it downloads. Features: avcodec, avformat, swscale, svt-av1 and
# none of gpl / nonfree / version3, so FFmpeg is LGPL 2.1 or later; the script fails otherwise.
# Needs Visual Studio 2022 (C++ workload) and git; about 30 minutes on a CI runner.
#
# Not verified locally (MDBridge's development machine is Linux): written for, and run by, CI.
param(
    [string]$Prefix = "build\ffmpeg\prefix",
    [string]$Work = "build\ffmpeg\work"
)
$ErrorActionPreference = "Stop"
$Here = $PSScriptRoot
$FFmpegDir = (Resolve-Path "$Here\..").Path

function Invoke-Checked {
    & $args[0] @($args | Select-Object -Skip 1)
    if ($LASTEXITCODE -ne 0) { throw "failed (exit $LASTEXITCODE): $($args -join ' ')" }
}

# pins: KEY=VALUE lines of versions.env
$pin = @{}
foreach ($line in Get-Content "$FFmpegDir\versions.env") {
    if ($line -match '^\s*([A-Z0-9_]+)=(.*)$') { $pin[$Matches[1]] = $Matches[2].Trim() }
}

New-Item -ItemType Directory -Force $Work | Out-Null
$Work = (Resolve-Path $Work).Path
if ($Work -match ' ') { throw "the work folder must not contain spaces (FFmpeg's build): $Work" }
New-Item -ItemType Directory -Force $Prefix | Out-Null
$Prefix = (Resolve-Path $Prefix).Path

Write-Output "== vcpkg $($pin.VCPKG_TAG)"
$Vcpkg = "$Work\vcpkg"
if (-not (Test-Path "$Vcpkg\.git")) {
    Invoke-Checked git -c advice.detachedHead=false clone --quiet --depth 1 --branch $pin.VCPKG_TAG $pin.VCPKG_GIT $Vcpkg
}
$head = (git -C $Vcpkg rev-parse HEAD).Trim()
if ($head -ne $pin.VCPKG_COMMIT) { throw "vcpkg is at $head, expected $($pin.VCPKG_COMMIT) ($($pin.VCPKG_TAG))" }
$ffPort = Get-Content "$Vcpkg\ports\ffmpeg\vcpkg.json" -Raw | ConvertFrom-Json
$svtPort = Get-Content "$Vcpkg\ports\svt-av1\vcpkg.json" -Raw | ConvertFrom-Json
if ($ffPort.version -ne $pin.FFMPEG_VERSION) { throw "vcpkg's ffmpeg port is $($ffPort.version), expected $($pin.FFMPEG_VERSION)" }
if ($svtPort.'version-semver' -ne $pin.SVTAV1_VERSION) {
    throw "vcpkg's svt-av1 port is $($svtPort.'version-semver'), expected $($pin.SVTAV1_VERSION)"
}
if (-not (Test-Path "$Vcpkg\vcpkg.exe")) {
    Invoke-Checked "$Vcpkg\bootstrap-vcpkg.bat" -disableMetrics
}

Write-Output "== FFmpeg $($pin.FFMPEG_VERSION) and SVT-AV1 $($pin.SVTAV1_VERSION) (x64-windows-release: DLLs, release only)"
$Triplet = "x64-windows-release"
$Installed = "$Work\installed"
$env:VCPKG_DISABLE_METRICS = "1"
Invoke-Checked "$Vcpkg\vcpkg.exe" install --classic --triplet $Triplet "--x-install-root=$Installed" `
    "ffmpeg[core,avcodec,avformat,swscale,svt-av1]"
$Root = "$Installed\$Triplet"

Write-Output "== check (LGPL only)"
$log = Get-ChildItem "$Vcpkg\buildtrees\ffmpeg" -Filter "build-$Triplet-rel-out.log" | Select-Object -First 1
$license = ""
if ($log) {
    $m = Select-String -Path $log.FullName -Pattern '^License: (.*)$' | Select-Object -First 1
    if ($m) { $license = $m.Matches[0].Groups[1].Value.Trim() }
}
if ($license -ne "LGPL version 2.1 or later") {
    throw "FFmpeg's configure reports the licence '$license'; MDBridge needs 'LGPL version 2.1 or later'"
}
$avcodec = Get-ChildItem "$Root\bin" -Filter "avcodec-*.dll" | Select-Object -First 1
if (-not $avcodec) { throw "no avcodec DLL in $Root\bin" }
# the configuration string avcodec_configuration() returns is a literal in the DLL
$text = [System.Text.Encoding]::ASCII.GetString([System.IO.File]::ReadAllBytes($avcodec.FullName))
$config = [regex]::Match($text, '--prefix=[\x20-\x7e]*').Value
if (-not $config) { throw "no configuration string in $($avcodec.Name)" }
foreach ($flag in "--enable-gpl", "--enable-nonfree", "--enable-version3") {
    if (" $config " -match " $([regex]::Escape($flag)) ") { throw "FFmpeg was configured with $flag`: $config" }
}
if ($config -notmatch '--enable-libsvtav1') { throw "FFmpeg was configured without libsvtav1: $config" }

Write-Output "== prefix $Prefix"
foreach ($dir in "include", "lib", "bin", "share\mdbridge-ffmpeg") {
    if (Test-Path "$Prefix\$dir") { Remove-Item -Recurse -Force "$Prefix\$dir" }
    New-Item -ItemType Directory -Force "$Prefix\$dir" | Out-Null
}
foreach ($lib in "libavcodec", "libavformat", "libavutil", "libswscale") {
    Copy-Item -Recurse "$Root\include\$lib" "$Prefix\include\"
}
foreach ($lib in "avcodec", "avformat", "avutil", "swscale") {
    Copy-Item "$Root\lib\$lib.lib" "$Prefix\lib\"
}
Copy-Item "$Root\bin\*.dll" "$Prefix\bin\"
$info = "$Prefix\share\mdbridge-ffmpeg"
Copy-Item "$Root\share\ffmpeg\copyright" "$info\FFmpeg-copyright.txt"
Copy-Item "$Root\share\svt-av1\copyright" "$info\SVT-AV1-copyright.txt"
@(
    "MDBridge video libraries, built by packaging/ffmpeg/windows/build_ffmpeg.ps1"
    "FFmpeg $($pin.FFMPEG_VERSION) (vcpkg port ffmpeg $($ffPort.version)#$($ffPort.'port-version'), source tag n$($pin.FFMPEG_VERSION))"
    "SVT-AV1 $($pin.SVTAV1_VERSION) (vcpkg port svt-av1 $($svtPort.'version-semver'))"
    "vcpkg $($pin.VCPKG_TAG) ($($pin.VCPKG_GIT), commit $($pin.VCPKG_COMMIT)), triplet $Triplet"
    "FFmpeg licence: $license"
    "FFmpeg configuration: $config"
    "DLLs: $((Get-ChildItem "$Prefix\bin" -Filter *.dll | ForEach-Object Name) -join ', ')"
) | Set-Content -Encoding utf8 "$info\BUILDINFO.txt"

Write-Output "== corresponding source of the FFmpeg DLLs"
$stage = "$Work\source-stage\ffmpeg-$($pin.FFMPEG_VERSION)-vcpkg"
if (Test-Path "$Work\source-stage") { Remove-Item -Recurse -Force "$Work\source-stage" }
New-Item -ItemType Directory -Force "$stage\vcpkg-port-ffmpeg", "$stage\MDBRIDGE-BUILD" | Out-Null
$archive = Get-ChildItem "$Vcpkg\downloads" -Filter "ffmpeg-ffmpeg-n$($pin.FFMPEG_VERSION)*.tar.gz" | Select-Object -First 1
if (-not $archive) { throw "FFmpeg's source archive is not in $Vcpkg\downloads" }
Copy-Item $archive.FullName "$stage\"
Copy-Item "$Vcpkg\ports\ffmpeg\*" "$stage\vcpkg-port-ffmpeg\"
Copy-Item "$FFmpegDir\versions.env", "$Here\build_ffmpeg.ps1", "$info\BUILDINFO.txt" "$stage\MDBRIDGE-BUILD\"
@(
    "Corresponding source of the FFmpeg DLLs bundled with MDBridge for Windows."
    "$($archive.Name): FFmpeg's source as vcpkg $($pin.VCPKG_TAG) downloaded it (tag n$($pin.FFMPEG_VERSION))."
    "vcpkg-port-ffmpeg: vcpkg's build recipe and the patches it applies (portfile.cmake lists them)."
    "MDBRIDGE-BUILD: the script that ran vcpkg, the pinned versions and the resulting configuration."
    "Rebuild: clone vcpkg at commit $($pin.VCPKG_COMMIT) and run build_ffmpeg.ps1 (needs Visual Studio 2022)."
) | Set-Content -Encoding utf8 "$stage\README.txt"
$zip = Join-Path (Split-Path $Prefix -Parent) "ffmpeg-$($pin.FFMPEG_VERSION)-vcpkg-source.zip"
Compress-Archive -Path $stage -DestinationPath $zip -Force
Remove-Item -Recurse -Force "$Work\source-stage"

Write-Output "== done"
Get-Content "$info\BUILDINFO.txt" | ForEach-Object { "   $_" }
Get-Item $zip | ForEach-Object { "   source: $($_.Name) ($([math]::Round($_.Length / 1MB, 1)) MB)" }

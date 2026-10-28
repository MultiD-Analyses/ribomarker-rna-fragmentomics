<#
ribomarker_windows.ps1: the RiboMarker pipeline in Docker, on the samples of a sample sheet.
Started by ribomarker.cmd, so that PowerShell's script setting (execution policy) does not need changing:

  ribomarker SHEET OUTPUT [options]    run the pipeline (SHEET: your sample sheet, a CSV file in the
                                       folder of the FASTQ files; OUTPUT: results)
  ribomarker --demo OUTPUT             test the installation (a few minutes)
  ribomarker --help                    all options

Options: --threads N (default 8), --species NAME, --run-id NAME, --config FILE,
--references FOLDER, --tiers FILE, --rna-classes FILE, --unlock. Folders and files are your own paths; this script mounts them into
the container (the sheet's folder and the files read-only) and writes the docker run command.
Image: $env:RIBOMARKER_IMAGE (default dockermultid/ribomarker-rna-fragmentomics:0.3.0). Same behaviour as the
bash script ribomarker (Linux, macOS, WSL).
#>
$ErrorActionPreference = 'Stop'
$Image = if ($env:RIBOMARKER_IMAGE) { $env:RIBOMARKER_IMAGE } else { 'dockermultid/ribomarker-rna-fragmentomics:0.3.0' }
$Docker = if ($env:RIBOMARKER_DOCKER) { $env:RIBOMARKER_DOCKER } else { 'docker' }

function Fail($message) {
    [Console]::Error.WriteLine("ribomarker: $message")
    exit 2
}
function FolderPath($path) { (Resolve-Path -LiteralPath $path).Path }

$dockerArgs = @('run', '--rm', '--init', '--network', 'none')
if (-not [Console]::IsInputRedirected -and -not [Console]::IsOutputRedirected) {
    $dockerArgs += '-it'            # a terminal: live output, and Ctrl-C stops the run cleanly
}
$runArgs = @()
$paths = @()
$demo = $false; $help = $false; $threads = '8'

$i = 0
while ($i -lt $args.Count) {
    $a = [string]$args[$i]
    if ($a -eq '-h' -or $a -eq '--help') { $help = $true; $i += 1 }
    elseif ($a -eq '--demo') { $demo = $true; $i += 1 }
    elseif ($a -eq '--unlock') { $runArgs += '--unlock'; $i += 1 }
    elseif ($a -in @('--threads', '--species', '--run-id')) {
        if ($i + 1 -ge $args.Count) { Fail "$a needs a value" }
        if ($a -eq '--threads') { $threads = [string]$args[$i + 1] } else { $runArgs += $a, [string]$args[$i + 1] }
        $i += 2
    }
    elseif ($a -in @('--config', '--tiers', '--rna-classes')) {     # a file: mounted read-only as /data/<option>/<name>
        if ($i + 1 -ge $args.Count) { Fail "$a needs a file" }
        $file = [string]$args[$i + 1]
        if (-not (Test-Path -LiteralPath $file -PathType Leaf)) { Fail "${a}: file not found: $file" }
        $name = Split-Path -Leaf $file
        $kind = $a.Substring(2)
        $dockerArgs += '-v', "$(FolderPath $file):/data/$kind/${name}:ro"
        $runArgs += $a, "/data/$kind/$name"
        $i += 2
    }
    elseif ($a -eq '--references') {
        if ($i + 1 -ge $args.Count) { Fail "$a needs a folder" }
        $folder = [string]$args[$i + 1]
        if (-not (Test-Path -LiteralPath $folder -PathType Container)) { Fail "${a}: folder not found: $folder" }
        $dockerArgs += '-v', "$(FolderPath $folder):/data/references:ro"
        $runArgs += '--references', '/data/references'
        $i += 2
    }
    elseif ($a.StartsWith('-')) { Fail "unknown option $a (ribomarker --help lists the options)" }
    else { $paths += $a; $i += 1 }
}

if (-not (Get-Command $Docker -ErrorAction SilentlyContinue)) {
    Fail 'Docker is not installed; see https://docs.docker.com/get-docker/'
}
# docker info can print harmless warnings (e.g. "No blkio throttle" on WSL2). Windows PowerShell 5.1
# turns redirected warnings of a program into errors, so errors must not stop the script during this check.
$ErrorActionPreference = 'Continue'
& $Docker info *> $null
$dockerRunning = ($LASTEXITCODE -eq 0)
$ErrorActionPreference = 'Stop'
if (-not $dockerRunning) { Fail 'Docker is not running; start Docker Desktop and try again' }

if ($help) {
    Write-Output 'ribomarker SHEET OUTPUT [options]   (SHEET: your sample sheet, in the folder of the FASTQ files; OUTPUT: the results folder; default 8 threads)'
    Write-Output 'ribomarker --demo OUTPUT            (test the installation)'
    Write-Output "The options are those of the pipeline's run command, with your own paths for files and folders:"
    Write-Output ''
    & $Docker run --rm $Image run --help
    exit $LASTEXITCODE
}

if ($demo) {
    if ($paths.Count -ne 1) { Fail 'usage: ribomarker --demo OUTPUT' }
    New-Item -ItemType Directory -Force -Path $paths[0] | Out-Null
    $dockerArgs += '-v', "$(FolderPath $paths[0]):/data/output"
    & $Docker @dockerArgs $Image run --demo /data/output --threads $threads
    exit $LASTEXITCODE
}

if ($paths.Count -ne 2) { Fail 'usage: ribomarker SHEET OUTPUT [options]  (ribomarker --help lists the options)' }
if (-not (Test-Path -LiteralPath $paths[0] -PathType Leaf)) { Fail "sample sheet not found: $($paths[0])" }
$sheet = FolderPath $paths[0]
$input_ = Split-Path -Parent $sheet                 # the sheet's folder, with the FASTQ files
$existed = Test-Path -LiteralPath $paths[1] -PathType Container
New-Item -ItemType Directory -Force -Path $paths[1] | Out-Null   # created here, so that Docker does not create it
$output = FolderPath $paths[1]
if (($output.TrimEnd('\') + '\').StartsWith($input_.TrimEnd('\') + '\', [StringComparison]::OrdinalIgnoreCase)) {
    if (-not $existed) { Remove-Item -LiteralPath $paths[1] }     # leave nothing behind in the input folder
    Fail 'the output folder must not be the folder of the sample sheet or lie inside it'
}
$dockerArgs += '-v', "${input_}:/data/input:ro", '-v', "${output}:/data/output"
& $Docker @dockerArgs $Image run "/data/input/$(Split-Path -Leaf $sheet)" /data/output --threads $threads @runArgs
exit $LASTEXITCODE

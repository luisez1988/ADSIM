<#
.SYNOPSIS
    Run a long command with the system kept awake, then let it sleep normally again.

.DESCRIPTION
    This machine sleeps after 5 minutes idle on AC and 3 minutes on battery, which stops
    a long ADSIM run within minutes of you walking away. A CPU-bound solver does NOT
    count as activity: Windows measures idleness by user input, not load.

    This wrapper asks Windows to hold the system awake for exactly as long as the child
    process runs, using SetThreadExecutionState(ES_CONTINUOUS | ES_SYSTEM_REQUIRED). It
    changes no setting and leaves nothing behind - the request dies with the process, so
    an interrupted or crashed run cannot leave the machine unable to sleep.

    The DISPLAY is still allowed to switch off; only the system is held awake. Pass
    -KeepDisplayOn if you want the screen up too (it costs battery and burns the panel,
    so it is off by default).

.PARAMETER Command
    Executable to run. Defaults to julia.

.PARAMETER Arguments
    Arguments for it. Everything after -- is passed through untouched.

.PARAMETER KeepDisplayOn
    Also keep the display awake.

.EXAMPLE
    # From the repo root - run the 2 mm case on 8 threads
    .\run_nosleep.ps1 -WorkingDirectory src -Arguments '-t8','kernel.jl','2mm_size'

.EXAMPLE
    # Anything else, not just ADSIM
    .\run_nosleep.ps1 -Command cmd -Arguments '/c','timeout','120'

.NOTES
    CLOSING THE LID STILL SLEEPS THE MACHINE. No process can override the lid switch.
    Either leave the lid open, or change the lid action first:
        powercfg /setacvalueindex SCHEME_CURRENT SUB_BUTTONS LIDACTION 0
        powercfg /setactive SCHEME_CURRENT
    (0 = do nothing; restore with 1 = sleep.)

    Stay on AC. On battery Windows will still sleep at a critical charge level, and the
    3-minute DC idle timeout is what this wrapper suppresses - not the battery itself.
#>

[CmdletBinding()]
param(
    [string]   $Command = 'julia',
    [string[]] $Arguments = @(),
    [string]   $WorkingDirectory = $PWD.Path,
    [switch]   $KeepDisplayOn
)

$ErrorActionPreference = 'Stop'

$signature = @'
[DllImport("kernel32.dll", SetLastError = true)]
public static extern uint SetThreadExecutionState(uint esFlags);
'@
$Power = Add-Type -MemberDefinition $signature -Name NoSleep -Namespace ADSIM -PassThru

$ES_CONTINUOUS       = [uint32]'0x80000000'
$ES_SYSTEM_REQUIRED  = [uint32]'0x00000001'
$ES_DISPLAY_REQUIRED = [uint32]'0x00000002'

$flags = $ES_CONTINUOUS -bor $ES_SYSTEM_REQUIRED
if ($KeepDisplayOn) { $flags = $flags -bor $ES_DISPLAY_REQUIRED }

# Assert once before starting so there is no window where the machine could drop off
# between here and the first heartbeat.
if ($Power::SetThreadExecutionState($flags) -eq 0) {
    throw 'SetThreadExecutionState failed; the system was NOT held awake. Aborting rather than starting a run that will be cut short.'
}

$started = Get-Date
Write-Host "Keeping system awake (display $(if ($KeepDisplayOn) { 'ON' } else { 'may sleep' }))." -ForegroundColor Cyan
Write-Host "Running: $Command $($Arguments -join ' ')" -ForegroundColor Cyan
Write-Host "  in: $WorkingDirectory"
Write-Host "  NOTE: closing the lid will still sleep this machine." -ForegroundColor Yellow

$proc = $null
try {
    $proc = Start-Process -FilePath $Command -ArgumentList $Arguments `
                          -WorkingDirectory $WorkingDirectory -NoNewWindow -PassThru

    # Re-assert on a heartbeat rather than trusting a single call. The execution state is
    # per-thread, and PowerShell does not promise that this script body stays on one
    # thread for hours; re-asserting makes that irrelevant. It also gives a visible sign
    # the wrapper is still alive during a run that prints nothing for long stretches.
    # Re-assert every minute (well inside the 3-minute battery idle timeout) but only
    # report every tenth, so a multi-hour run does not bury the solver's own output.
    $beat = 0
    while (-not $proc.WaitForExit(60000)) {
        [void]$Power::SetThreadExecutionState($flags)
        $beat++
        if ($beat % 10 -eq 0) {
            $elapsed = (Get-Date) - $started
            Write-Host ("  [awake] {0:hh\:mm\:ss} elapsed, pid {1} running" -f $elapsed, $proc.Id) -ForegroundColor DarkGray
        }
    }
}
finally {
    # Release unconditionally: a crash, a Ctrl+C or a throw must not leave the machine
    # pinned awake. ES_CONTINUOUS with no other flag clears the standing request.
    [void]$Power::SetThreadExecutionState($ES_CONTINUOUS)
    Write-Host 'Released the wake request; normal sleep behaviour restored.' -ForegroundColor Cyan
}

$elapsed = (Get-Date) - $started
Write-Host ("Finished in {0:hh\:mm\:ss} with exit code {1}." -f $elapsed, $proc.ExitCode) `
           -ForegroundColor $(if ($proc.ExitCode -eq 0) { 'Green' } else { 'Red' })
exit $proc.ExitCode

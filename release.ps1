# release.ps1 - Push to GitHub (with Release) and PyPI in one shot
# Usage: .\release.ps1 "your commit message"

param(
    [Parameter(Mandatory=$true)]
    [string]$Message
)

$ErrorActionPreference = "Stop"

$version = (Select-String -Path "pyproject.toml" -Pattern 'version = "(.+)"').Matches.Groups[1].Value
# Preflight: all version files must agree before anything ships (v0.7.2 lesson)
$initV = (Select-String -Path "src/mtft/__init__.py" -Pattern '__version__ = "(.+)"').Matches.Groups[1].Value
$cffV  = (Select-String -Path "CITATION.cff" -Pattern '^version:\s*(\S+)').Matches.Groups[1].Value
if (($initV -ne $version) -or ($cffV -ne $version)) {
    Write-Host "VERSION MISMATCH — fix before releasing:" -ForegroundColor Red
    Write-Host "  pyproject.toml : $version"  -ForegroundColor Red
    Write-Host "  __init__.py    : $initV"    -ForegroundColor Red
    Write-Host "  CITATION.cff   : $cffV"     -ForegroundColor Red
    exit 1
}
Write-Host "Releasing mtft v$version" -ForegroundColor Cyan
Write-Host "Message: $Message" -ForegroundColor Gray
Write-Host ""

# Step 1: Git commit and push
Write-Host "-- Git -------------------------------------------------" -ForegroundColor Yellow
git add -A
git commit -m "$Message"
git push origin main
Write-Host "Git push: done" -ForegroundColor Green
Write-Host ""

# Step 2: Build Python package
Write-Host "-- Build -----------------------------------------------" -ForegroundColor Yellow
Remove-Item -Recurse -Force dist, build -ErrorAction SilentlyContinue
Remove-Item -Recurse -Force src\mtft.egg-info -ErrorAction SilentlyContinue
py -m build
Write-Host "Build: done" -ForegroundColor Green
Write-Host ""

# Step 3: Create GitHub Release (includes tag + release page)
Write-Host "-- GitHub Release --------------------------------------" -ForegroundColor Yellow
gh release create "v$version" dist\* --title "v$version" --notes "$Message" --latest
Write-Host "GitHub Release: done" -ForegroundColor Green
Write-Host ""

# Step 4: Verify wheel
Write-Host "-- Verify ----------------------------------------------" -ForegroundColor Yellow
py -m twine check dist\*
Write-Host ""

# Step 5: PyPI publish is handled by GitHub Actions (trusted publishing)
Write-Host "-- PyPI ------------------------------------------------" -ForegroundColor Yellow
Write-Host "GitHub Actions is re-running the full test suite on a clean machine" -ForegroundColor Gray
Write-Host "and publishing to PyPI. Nothing is published if any test fails." -ForegroundColor Gray
Write-Host "Watch: https://github.com/Kaizoku-Ronin/mtft/actions/workflows/publish.yml" -ForegroundColor Gray
Write-Host ""

Write-Host "========================================================" -ForegroundColor Cyan
Write-Host "  mtft v$version released" -ForegroundColor Green
Write-Host "  PyPI:   https://pypi.org/project/mtft/$version/  (live when the workflow finishes)" -ForegroundColor Gray
Write-Host "  GitHub: https://github.com/Kaizoku-Ronin/mtft/releases/tag/v$version" -ForegroundColor Gray
Write-Host "========================================================" -ForegroundColor Cyan

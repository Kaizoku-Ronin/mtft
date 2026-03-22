# release.ps1 - Push to GitHub and PyPI in one shot
# Usage: .\release.ps1 "your commit message"

param(
    [Parameter(Mandatory=$true)]
    [string]$Message
)

$ErrorActionPreference = "Stop"

$version = (Select-String -Path "pyproject.toml" -Pattern 'version = "(.+)"').Matches.Groups[1].Value
Write-Host "Releasing mtft v$version" -ForegroundColor Cyan
Write-Host "Message: $Message" -ForegroundColor Gray
Write-Host ""

# Step 1: Git
Write-Host "-- Git -------------------------------------------------" -ForegroundColor Yellow
git add -A
git commit -m "$Message"
git tag -a "v$version" -m "Release v$version"
git push origin main
git push origin "v$version"
Write-Host "Git: done" -ForegroundColor Green
Write-Host ""

# Step 2: Build
Write-Host "-- Build -----------------------------------------------" -ForegroundColor Yellow
Remove-Item -Recurse -Force dist, build -ErrorAction SilentlyContinue
Remove-Item -Recurse -Force src\mtft.egg-info -ErrorAction SilentlyContinue
py -m build
Write-Host "Build: done" -ForegroundColor Green
Write-Host ""

# Step 3: Verify
Write-Host "-- Verify ----------------------------------------------" -ForegroundColor Yellow
py -m twine check dist\*
Write-Host ""

# Step 4: Upload to PyPI
Write-Host "-- Upload to PyPI --------------------------------------" -ForegroundColor Yellow
py -m twine upload dist\*
Write-Host ""

Write-Host "========================================================" -ForegroundColor Cyan
Write-Host "  mtft v$version is live on GitHub + PyPI" -ForegroundColor Green
Write-Host "  https://pypi.org/project/mtft/$version/" -ForegroundColor Gray
Write-Host "========================================================" -ForegroundColor Cyan

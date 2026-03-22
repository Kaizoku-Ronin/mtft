# release.ps1 - Push to GitHub (with Release) and PyPI in one shot
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

# Step 5: Upload to PyPI
Write-Host "-- Upload to PyPI --------------------------------------" -ForegroundColor Yellow
py -m twine upload dist\*
Write-Host ""

Write-Host "========================================================" -ForegroundColor Cyan
Write-Host "  mtft v$version is live everywhere" -ForegroundColor Green
Write-Host "  PyPI:   https://pypi.org/project/mtft/$version/" -ForegroundColor Gray
Write-Host "  GitHub: https://github.com/Kaizoku-Ronin/mtft/releases/tag/v$version" -ForegroundColor Gray
Write-Host "========================================================" -ForegroundColor Cyan

# Trusted Publishing Setup — one-time, ~2 minutes

Trusted Publishing links your GitHub repo directly to your PyPI project.
After setup, creating a GitHub Release automatically tests, builds, and
publishes to PyPI. No API tokens anywhere — GitHub proves its identity
to PyPI cryptographically (OpenID Connect) on every release.

## One-time setup on PyPI

1. Log in at https://pypi.org as **Kaizoku_Ronin**
2. Go to **Your projects → mtft → Settings** (Manage project)
3. Scroll to **Publishing** and click **Add a new publisher**
4. Choose the **GitHub** tab and enter exactly:

   | Field | Value |
   |---|---|
   | Owner | `Kaizoku-Ronin` |
   | Repository name | `mtft` |
   | Workflow name | `publish.yml` |
   | Environment name | `pypi` |

5. Click **Add**

That's it. Nothing to configure on the GitHub side — the workflow file
(`.github/workflows/publish.yml`) is already in the repo, and the
`environment: pypi` in it matches what you entered above.

## Releasing from now on

Your existing habit stays the same:

```powershell
.\release.ps1 "v0.6.1 — audit coalescence release"
```

The script commits, pushes, and creates the GitHub Release. The Release
triggers the workflow, which:

1. runs all 235 tests on a clean machine (publish is BLOCKED if any fail)
2. builds the sdist + wheel
3. validates with twine check
4. publishes to PyPI

Watch progress at: https://github.com/Kaizoku-Ronin/mtft/actions

The `py -m twine upload` step in release.ps1 is now redundant — you can
delete Step 5 from the script (and your PyPI API token, if you like).
If you leave Step 5 in, the worst case is a harmless "File already
exists" error when twine races the workflow.

## Why this is better than the token

- Your PyPI token no longer needs to exist on your laptop at all — one
  less credential to protect (aligned with the MonsterShield practice
  of not leaving key material lying around).
- Every published artifact is built on a clean GitHub runner from the
  exact tagged commit — reproducible, and the test gate means a broken
  build can never reach PyPI.
- PyPI marks these releases with a verified provenance badge.

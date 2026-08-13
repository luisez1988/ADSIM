# Quick Start Guide: Your First GitHub Release

## You're Ready! 🎉

All version management and release automation is now set up. Here's what to do next:

---

## Option 1: Create Release Now (Recommended)

If you want to create the v0.1.0 release immediately:

### Step 1: Commit Everything

```powershell
# Navigate to project root
cd c:\Users\zamcr\Dcuments\GitHub\ADSIM

# Check what will be committed
git status

# Add all new files
git add .

# Commit with descriptive message
git commit -m "Add version management and release automation system

- Add centralized VERSION file (0.1.0)
- Add automated version update script
- Add GitHub Actions workflows for CI/CD
- Add comprehensive CHANGELOG.md
- Update README with installation and release info
- Integrate version display in application logs and outputs
- Add release documentation and guides"
```

### Step 2: Create the Tag

```powershell
# Create annotated tag for v0.1.0
git tag -a v0.1.0 -m "First alpha release of ADSIM v0.1.0

Core Features:
- Multi-gas transport simulation
- Fully explicit finite element solver
- GiD problem type integration
- VTK output for ParaView
- Example test cases
- Windows standalone executable

This is an alpha release for testing and feedback."
```

### Step 3: Push to GitHub

```powershell
# Push the commits
git push origin main

# Push the tag (this triggers the automated build!)
git push origin v0.1.0
```

### Step 4: Monitor Progress

1. **Open your browser and go to:**
   ```
   https://github.com/luisez1988/ADSIM/actions
   ```

2. **You should see a "Release Build" workflow running**
   - Click on it to see real-time progress
   - It will take about 10-15 minutes to complete
   - Wait for the green checkmark ✅

3. **Once complete, visit:**
   ```
   https://github.com/luisez1988/ADSIM/releases
   ```

4. **Your release is live!** 
   - Download the ZIP file to test
   - Share the link with users

---

## Option 2: Test Everything First (Safer)

If you want to test locally before releasing:

### Build Executable Locally

```powershell
# Navigate to project root
cd c:\Users\zamcr\Documents\GitHub\ADSIM

# One-time setup: resolve the package environment and install the build tool
julia --project=. -e "using Pkg; Pkg.instantiate()"
julia -e "using Pkg; Pkg.add(\"PackageCompiler\")"

# Build the executable
julia build_app.jl

# This will take 10-15 minutes...
```

The build requires, in the project root:

- `Project.toml` declaring the package name, uuid and dependencies, plus the
  `Manifest.toml` produced by `Pkg.instantiate()`. `create_app` builds a package,
  not a loose folder of scripts.
- `julia_main()` in `src/ADSIM.jl`. This is the entry point PackageCompiler calls
  in the compiled app; `main()` alone is not enough.
- `PackageCompiler` installed in the environment that runs `build_app.jl`.

`buildscripts/run_cli.jl` traces compilation by running a real analysis. It copies
its case into a temporary workspace, so a build never writes VTK files, probe
series or checkpoints into `src/output/`.

### Test the Executable

```powershell
# Test version display
.\ADSIM_app\bin\ADSIM.exe --version
# Should output: ADSIM v0.1.0

# Test with example problem
cd src
..\ADSIM_app\bin\ADSIM.exe Pani_etal

# Check that it runs successfully
```

### Verify Output Files

```powershell
# Check log file contains version
Select-String -Path "output\Pani_etal.log" -Pattern "Version:"
# Should show: Version: 0.1.0

# Check VTK files exist
Get-ChildItem output\*.vtk | Select-Object Name

# Check VTK header has version
Get-Content output\Pani_etal_000000.vtk -Head 5
# Should show: ADSIM v0.1.0 - Time Step 0, Time = 0.0
```

### If All Tests Pass

Then proceed with **Option 1** above to create the release!

---

## What Happens When You Push the Tag?

GitHub Actions will automatically:

1. ✅ **Checkout your code** - Gets the latest version from the tag
2. ✅ **Setup Julia 1.8** - Installs Julia environment
3. ✅ **Install dependencies** - Runs `Pkg.instantiate()`
4. ✅ **Build executable** - Uses PackageCompiler to create ADSIM.exe
5. ✅ **Package as ZIP** - Creates `ADSIM-v0.1.0-windows-x64.zip`
6. ✅ **Create release** - Makes a new GitHub release
7. ✅ **Upload artifact** - Attaches ZIP to the release
8. ✅ **Extract release notes** - From CHANGELOG.md

You don't need to do anything manually - it's all automated!

---

## After the Release

### Share Your Release

Once the release is published, you can share:

**Direct download link:**
```
https://github.com/luisez1988/ADSIM/releases/download/v0.1.0/ADSIM-v0.1.0-windows-x64.zip
```

**Release page:**
```
https://github.com/luisez1988/ADSIM/releases/tag/v0.1.0
```

### Update Your README Badge

The version badge will automatically show v0.1.0:
```markdown
[![Version](https://img.shields.io/badge/version-0.1.0-blue.svg)](https://github.com/luisez1988/ADSIM/releases)
```

---

## Future Releases

For version 0.2.0 or later:

```powershell
# 1. Update CHANGELOG.md with new changes

# 2. Update version
julia scripts/update_version.jl 0.2.0

# 3. Commit and tag
git add .
git commit -m "Bump version to v0.2.0"
git tag -a v0.2.0 -m "Release v0.2.0 - New features"
git push origin main
git push origin v0.2.0

# 4. Watch the magic happen at:
# https://github.com/luisez1988/ADSIM/actions
```

---

## Need Help?

Refer to these guides:

- **Complete guide:** `.github/HOW_TO_RELEASE.md`
- **Release checklist:** `.github/RELEASE_TEMPLATE.md`
- **Version management:** `VERSION_MANAGEMENT.md`
- **Changelog:** `CHANGELOG.md`

---

## Summary of What Was Created

### New Files (14 total)

1. `VERSION` - Single source of truth (0.1.0)
2. `src/version.jl` - Version reader module
3. `scripts/update_version.jl` - Automated version updater
4. `CHANGELOG.md` - Version history
5. `.github/workflows/release.yml` - Release automation
6. `.github/workflows/test.yml` - CI testing
7. `.github/RELEASE_TEMPLATE.md` - Release checklist
8. `.github/HOW_TO_RELEASE.md` - Complete release guide
9. `VERSION_MANAGEMENT.md` - System documentation
10. `QUICKSTART_RELEASE.md` - This file!

### Modified Files (6 total)

1. `README.md` - Added badges and installation instructions
2. `src/ADSIM.jl` - Integrated version module
3. `src/kernel.jl` - Added --version flag and version display
4. `src/write_vtk.jl` - Added version to VTK headers
5. `Project.toml` - Version updated to 0.1.0
6. All GiD files - Version updated to 0.1.0

### Version Appears In:

- ✅ Log file banner
- ✅ Command line (`--version` flag)
- ✅ VTK file headers
- ✅ GiD splash screen
- ✅ Generated TOML/mesh files
- ✅ GitHub release page

---

## You're All Set! 🚀

Choose your path:
- **Ready to go?** → Follow **Option 1** (commit, tag, push)
- **Want to test first?** → Follow **Option 2** (build locally, then release)

Either way, you now have a professional version management and release system!

Good luck with your first release! 🎉

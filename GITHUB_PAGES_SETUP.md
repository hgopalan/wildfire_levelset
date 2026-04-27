# GitHub Pages Setup Instructions

## Issue
The documentation link https://hgopalan.github.io/wildfire_levelset/ is returning a 404 error because GitHub Pages deployment is failing.

## Root Cause
The GitHub Actions workflow is building the documentation successfully, but the deployment step fails with:
```
Error: Failed to create deployment (status: 404)
Ensure GitHub Pages has been enabled: https://github.com/hgopalan/wildfire_levelset/settings/pages
```

This indicates that GitHub Pages needs to be configured to use **GitHub Actions** as the deployment source.

## Solution

### Step 1: Configure GitHub Pages Source

1. Go to your repository settings: https://github.com/hgopalan/wildfire_levelset/settings/pages

2. Under "Build and deployment":
   - **Source**: Select **"GitHub Actions"** (not "Deploy from a branch")
   
   This is crucial! The workflow uses the modern `actions/deploy-pages@v4` action which requires GitHub Actions as the source.

### Step 2: Trigger a New Build (Optional)

Once the above is configured, you can either:

**Option A: Wait for the next push to main**
- The workflow will run automatically on the next commit to `main`

**Option B: Manually trigger the workflow**
1. Go to Actions tab: https://github.com/hgopalan/wildfire_levelset/actions
2. Select "Build and Deploy Documentation" workflow
3. Click "Run workflow" button
4. Select `main` branch
5. Click "Run workflow"

**Option C: Make an empty commit to trigger the workflow**
```bash
git checkout main
git commit --allow-empty -m "Trigger documentation rebuild"
git push origin main
```

### Step 3: Verify Deployment

After the workflow completes successfully:
- Check the Actions tab to ensure both "build" and "deploy" jobs completed successfully
- Visit https://hgopalan.github.io/wildfire_levelset/ (may take a few minutes to propagate)

## Current Status

✅ **Documentation source files**: All `.rst` files are present in `/docs`
✅ **Sphinx configuration**: Properly configured in `docs/conf.py`  
✅ **GitHub Actions workflow**: Workflow file `.github/workflows/docs.yml` is correct and active
✅ **Documentation builds successfully**: HTML pages generate without errors
❌ **Deployment**: Fails due to GitHub Pages not configured for Actions deployment

## Verification

I have successfully built the HTML documentation locally in this environment:
- HTML files are in `docs/_build/html/`
- No errors or warnings (except a harmless `_static` warning which is normal)
- All pages generated correctly:
  - index.html
  - overview.html
  - mathematical_models.html
  - code_structure.html
  - building.html
  - usage.html
  - api_reference.html

## Workflow Details

The workflow (`.github/workflows/docs.yml`):
1. ✅ Triggers on push to main/master
2. ✅ Has correct permissions (pages: write, id-token: write)
3. ✅ Builds documentation with Sphinx
4. ✅ Uploads artifact
5. ❌ Deploys to GitHub Pages (fails due to configuration)

The last workflow run (ID: 25022335971):
- Build job: ✅ Success
- Deploy job: ❌ Failed with "Not Found" error

## Next Steps

**For the repository owner:**
Configure GitHub Pages source to "GitHub Actions" in repository settings (see Step 1 above).

**Once configured:**
The documentation will be automatically built and deployed on every push to `main` branch.

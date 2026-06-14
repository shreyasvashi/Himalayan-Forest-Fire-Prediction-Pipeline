# Publishing This Project to GitHub

This guide explains how to organize the project files and push them to a new GitHub repository.

## 1. Final Folder Layout

Arrange your local project folder exactly like this before pushing:

```
himalayan-fire-prediction/
  data/
    .gitkeep
  output/
    .gitkeep
  fire_prediction_pipeline.py
  requirements.txt
  README.md
  .gitignore
  LICENSE
```

Notes:

- The `data/.gitkeep` and `output/.gitkeep` files are empty placeholder files. Git does not track empty directories, so these placeholders keep the `data` and `output` folders present in the repository even though their real contents (DEM files, NetCDF files, CSV files, and generated rasters) are excluded by `.gitignore`.
- Do not commit the actual DEM, NetCDF, or FIRMS CSV files, and do not commit generated GeoTIFF or JSON outputs. These files can be very large and may also be subject to redistribution restrictions from their original data providers. The `.gitignore` file already excludes these.
- Add a `LICENSE` file appropriate for your project. If you are not sure which license to use, the MIT license or the Apache 2.0 license are common choices for open source scientific tools.

## 2. Create the Repository on GitHub

1. Sign in to GitHub.
2. Click the plus icon in the top right corner and choose "New repository".
3. Enter a repository name, for example `himalayan-fire-prediction`.
4. Add a short description, for example "Cellular Automata fire spread prediction for the Himalayan region, calibrated with MCMC and driven by SRTM, ERA5-Land, and FIRMS data".
5. Choose "Public" or "Private" depending on your needs.
6. Do not initialize the repository with a README, .gitignore, or license at this stage, since you already have these files locally. If you prefer, you can initialize with these and merge afterward, but starting from an empty repository avoids merge conflicts.
7. Click "Create repository".

## 3. Initialize Git Locally and Push

Open a terminal in your project folder (the one containing `fire_prediction_pipeline.py`, `README.md`, and the other files) and run the following commands:

```
git init
git add .
git commit -m "Initial commit: Himalayan forest fire prediction pipeline"
git branch -M main
git remote add origin https://github.com/YOUR_USERNAME/himalayan-fire-prediction.git
git push -u origin main
```

Replace `YOUR_USERNAME` with your GitHub username, and replace `himalayan-fire-prediction` with the repository name you chose.

If you are prompted for credentials and you have two factor authentication enabled on your GitHub account, you will need to use a personal access token instead of your account password, or use SSH based authentication instead of HTTPS.

## 4. Using SSH Instead of HTTPS (Optional)

If you prefer SSH authentication:

1. Generate an SSH key pair if you do not already have one:

```
ssh-keygen -t ed25519 -C "your_email@example.com"
```

2. Add the generated public key (the contents of `~/.ssh/id_ed25519.pub`) to your GitHub account under Settings, then SSH and GPG keys.
3. Use the SSH remote URL instead of the HTTPS one:

```
git remote add origin git@github.com:YOUR_USERNAME/himalayan-fire-prediction.git
git push -u origin main
```

## 5. Verifying the Upload

After pushing, refresh the repository page on GitHub. You should see:

- `fire_prediction_pipeline.py`
- `requirements.txt`
- `README.md`
- `.gitignore`
- `LICENSE`
- The `data` and `output` folders, each containing only a `.gitkeep` file

The README will be rendered automatically on the repository's main page.

## 6. Recommended Repository Settings

- Add topics to the repository such as `wildfire`, `cellular-automata`, `mcmc`, `geospatial`, `remote-sensing`, and `himalaya` to make the project easier to discover.
- Enable "Issues" if you want to track bugs or feature requests.
- Consider adding a `CONTRIBUTING.md` file if you expect external contributions, describing coding style, how to run tests, and how to submit pull requests.
- Consider adding a GitHub Actions workflow under `.github/workflows/` to run linting or basic syntax checks automatically on each push. A minimal example workflow file is shown below.

## 7. Optional: Basic GitHub Actions Workflow

Create the file `.github/workflows/lint.yml` with the following content to automatically check that the code compiles without syntax errors on each push:

```yaml
name: Lint

on: [push, pull_request]

jobs:
  syntax-check:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - name: Set up Python
        uses: actions/setup-python@v5
        with:
          python-version: "3.11"
      - name: Compile check
        run: python -m py_compile fire_prediction_pipeline.py
```

This workflow does not install the full geospatial dependency stack, since some of those packages require system level GDAL installation and can slow down continuous integration runs significantly. It only verifies that the Python source file has no syntax errors. If you want a more thorough check, extend the workflow to install GDAL and the packages listed in `requirements.txt`, then run a small smoke test against sample data.

## 8. Keeping Large Data Files Out of Git History

If you accidentally commit a large DEM or NetCDF file, removing it from the working tree alone is not enough, since it remains in the git history. Use a tool such as `git filter-repo` or the BFG Repo-Cleaner to remove large files from history before pushing, or simply start a fresh repository if the project is new and the commit has not yet been pushed.

For ongoing large data management, consider using Git LFS (Large File Storage) for any sample data files you do want to version, by running:

```
git lfs install
git lfs track "*.tif"
git lfs track "*.nc"
git add .gitattributes
```

Then commit and push as usual. Note that Git LFS has storage and bandwidth limits on free GitHub accounts, so this is best suited for small sample datasets used for testing rather than full resolution regional data.

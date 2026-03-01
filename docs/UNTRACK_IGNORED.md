# Remove Already-Tracked Ignored Files from Git

These files were committed *before* .gitignore. Run these commands to stop tracking them (files stay on disk):

```powershell
cd "c:\Users\praveeeee\Desktop\qq protein\quantum_protein_v3-3"

# 1. Stop tracking venv/ (~29k pip files — this is what caused the huge push)
git rm -r --cached venv/

# 2. Stop tracking __pycache__/ (run each if the folder exists)
git rm -r --cached core/__pycache__/
git rm -r --cached core/disca/__pycache__/
git rm -r --cached core/yopo/__pycache__/
git rm -r --cached quantum/__pycache__/
git rm -r --cached quantum/hamiltonian/__pycache__/
git rm -r --cached quantum/mapping/__pycache__/
git rm -r --cached quantum/noise/__pycache__/
git rm -r --cached quantum/optimization/__pycache__/
git rm -r --cached quantum/vqe/__pycache__/
git rm -r --cached thermodynamics/__pycache__/
git rm -r --cached utils/__pycache__/

# 3. (Optional) Stop tracking PDB files — users can download via pipeline
git rm --cached data/raw/2I9M.pdb data/raw/2JOF.pdb

# 4. Commit and push
git add .gitignore
git commit -m "chore: remove venv, __pycache__, PDB from tracking; fix .gitignore"
git push origin main
```

After this, future `git add` will respect .gitignore and pushes will be small.

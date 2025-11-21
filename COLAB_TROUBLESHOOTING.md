# Google Colab Troubleshooting Guide

## Issue: ModuleNotFoundError for model_config

If you see this error:
```
ModuleNotFoundError: No module named 'model_config'
```

### Solution 1: Reload the Notebook (Recommended)

Colab may be caching an old version of the notebook. To get the latest version:

1. In Colab, go to **File** → **Revert to saved**
2. Or close the tab and reopen the notebook from GitHub:
   - Go to [Google Colab](https://colab.research.google.com/)
   - Click **File** → **Open notebook** → **GitHub** tab
   - Enter: `jkrantz159/XeNH`
   - Make sure branch is: `claude/google-colab-support-01CL1eMMoLYfXqHoiGmoNCdM`
   - Select: `XeNH_Simulation_Colab.ipynb`

3. Make sure to run cells in order, especially the **Verify Setup** cell

### Solution 2: Manual sys.path Fix

If reloading doesn't work, insert this cell **BEFORE** the simulation cell:

```python
# Manual fix for imports
import sys
import os

current_dir = os.getcwd()
print(f"Current directory: {current_dir}")

if current_dir not in sys.path:
    sys.path.insert(0, current_dir)
    print(f"✓ Added {current_dir} to sys.path")

# Verify files exist
import os
files = ['model_config.py', 'model_utils.py', 'parallel_xe_model.py', 'parallel_n_model.py']
for f in files:
    print(f"{'✓' if os.path.exists(f) else '✗'} {f}")

# Test import
try:
    from model_config import ModelConfig
    print("✓ Import successful!")
except ImportError as e:
    print(f"✗ Import failed: {e}")
    print("\nDirectory contents:")
    !ls -la
    print(f"\nPython path: {sys.path[:3]}")
```

### Solution 3: Check Clone Step

Make sure the clone cell output shows:
```
✓ Repository cloned from branch: claude/google-colab-support-01CL1eMMoLYfXqHoiGmoNCdM
```

And that you're in the correct directory:
```
/content/XeNH/python/src
```

If you see a different directory, the `%cd` command may have failed. Try running:
```python
import os
os.chdir('/content/XeNH/python/src')
print(f"Current directory: {os.getcwd()}")
!ls -la
```

### Solution 4: Complete Fresh Start

If nothing else works:

1. **Runtime** → **Disconnect and delete runtime**
2. **Runtime** → **Run all** (runs all cells from scratch)
3. Watch for any errors in the clone or verify steps

## Common Issues

### Issue: "Repository already exists" but files missing

The clone step detected an existing `XeNH` folder but it may be from a different branch. Fix:

```python
!rm -rf /content/XeNH
!git clone -b claude/google-colab-support-01CL1eMMoLYfXqHoiGmoNCdM https://github.com/jkrantz159/XeNH.git
%cd XeNH/python/src
!ls -la
```

### Issue: Wrong branch cloned

Check the branch with:
```python
%cd /content/XeNH
!git branch -a
!git log --oneline -5
```

Should show `claude/google-colab-support-01CL1eMMoLYfXqHoiGmoNCdM` and recent commits including "Add Python implementation".

## Verification Checklist

Before running the simulation, verify:

- [ ] Clone completed successfully
- [ ] Current directory is `/content/XeNH/python/src`
- [ ] All 4 required .py files exist (model_config.py, model_utils.py, parallel_xe_model.py, parallel_n_model.py)
- [ ] `sys.path` includes current directory
- [ ] Test import of `ModelConfig` succeeds

The **Verify Setup** cell (cell 6) should check all of these automatically.

## Still Having Issues?

If none of these solutions work:

1. Copy the output from the **Verify Setup** cell
2. Share it along with the exact error message
3. Check if the branch has the Python files:
   - Visit: https://github.com/jkrantz159/XeNH/tree/claude/google-colab-support-01CL1eMMoLYfXqHoiGmoNCdM/python/src
   - Verify you can see model_config.py and other files

## After Main Branch Merge

Once the Python implementation is merged into the main branch, the clone command will be simplified to:
```python
!git clone https://github.com/jkrantz159/XeNH.git
```

And these import issues should not occur.

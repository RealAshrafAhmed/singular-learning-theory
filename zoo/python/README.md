# Zoo Python Research Monorepo

## Quick Start

```bash
# Install all packages
uv sync --all-packages

# Run Jupyter Lab
uv run jupyter lab
```

## Common Commands

### Packages

```bash
# Install/sync all packages
uv sync --all-packages

# Add a dependency to a specific lib
uv add numpy --package sklearn_extensions

# Add a dev dependency to root
uv add pytest --dev

# Reinstall a specific package (after changes)
uv sync --all-packages --reinstall-package sklearn-extensions

# Force full reinstall (nuclear option)
rm uv.lock && uv sync --all-packages
```

### Jupyter

```bash
# Start Jupyter Lab
uv run jupyter lab

# Register kernel (if kernel issues)
uv run python -m ipykernel install --user --name=zoo-python
```

**Note:** Always restart the Jupyter kernel after installing/updating packages.

### Development

```bash
# Format code
uv run poe fmt

# Lint code
uv run poe lint

# Type check
uv run poe check

# Run tests
uv run poe test

# Run all checks
uv run poe all
```

### Creating New Packages

```bash
# Create a new library
uv init libs/newlib --lib

# Create a new project
uv init projects/newproject --package
```

Then run `uv sync --all-packages` to register it.

## Project Structure

```
├── libs/                   # Shared libraries
│   └── sklearn_extensions/ # Example: sklearn utilities
│       ├── pyproject.toml
│       └── src/sklearn_extensions/
├── projects/               # Research projects
│   ├── orca/
│   ├── luna/
│   └── tiny/
├── notebooks/              # Shared notebooks
├── pyproject.toml          # Root workspace config
└── uv.lock                 # Lockfile (committed)
```

## Troubleshooting

### Import not working in notebook

1. Restart the Jupyter kernel
2. If still broken: `uv sync --all-packages --reinstall-package <package-name>`
3. Check kernel is using correct venv:
   ```python
   import sys
   print(sys.executable)  # Should be .venv/bin/python
   ```

### Package not found after adding

Make sure it's in the workspace members in root `pyproject.toml`:

```toml
[tool.uv.workspace]
members = [
    "libs/*",
    "projects/orca",
    # ...
]
```

### Wrong Python version

```bash
uv python install 3.12
uv sync --all-packages
```
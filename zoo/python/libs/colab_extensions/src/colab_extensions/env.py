from pathlib import Path
import sys

def setup_env(repodir="/content/waterloo-slt-reading-group"):
    try:
        import google.colab
        IN_COLAB = True
    except ImportError:
        IN_COLAB = False

    if IN_COLAB:
        ROOT = Path(repodir)
    else:
        path = Path.cwd()
        while path != path.parent:
            if (path / "pyproject.toml").exists():
                ROOT = path
                break
            path = path.parent
        else:
            raise FileNotFoundError("Could not find project root")
    
    return {
        "IN_COLAB": IN_COLAB,
        "ROOT": ROOT,
        "DATA": ROOT / "data",
        "OUTPUTS": ROOT / "outputs",
    }
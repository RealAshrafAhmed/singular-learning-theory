import matplotlib.pyplot as plt


def configure_latex_plotting(font_size=12):
    """
    Configures Matplotlib to use the LaTeX backend for rendering text,
    ensuring high-quality mathematical symbols and formatting.

    NOTE: This requires the 'texlive-latex-base' and 'lmodern' packages
    to be installed on the system (usually via apt-get in Colab) AND
    the Colab runtime must have been restarted after installation.

    Args:
        font_size (int): The base font size to use for general text.
    """
    try:
        plt.rcParams.update(
            {
                "text.usetex": True,  # Enable the use of LaTeX for all text
                "font.family": "serif",  # Use a serif font (Latin Modern)
                "font.serif": ["Latin Modern Roman"],  # Specify the Latin Modern font
                "font.size": font_size,
                "axes.labelsize": font_size + 2,  # Axis labels slightly larger
                "xtick.labelsize": font_size,
                "ytick.labelsize": font_size,
                "legend.fontsize": font_size,
                "figure.titlesize": font_size + 4,
            }
        )
        print("Matplotlib successfully configured for LaTeX plotting.")
        print(f"Base font size set to: {font_size}")

    except Exception as e:
        print("ERROR: Failed to configure Matplotlib for LaTeX.")
        print(
            "This usually means the TeX packages are not installed or the runtime was not restarted."
        )
        print(f"Details: {e}")


# Example of configuration parameters you can set
DEFAULT_CONFIG = {
    "font_size": 12,
    # Add other default parameters here if needed
}

if __name__ == "__main__":
    configure_latex_plotting(**DEFAULT_CONFIG)

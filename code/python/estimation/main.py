import typer
from .mixnorm.ksigma import main as mixnorm_ksigma

app = typer.Typer(
    help="My powerful command-line interface application.",
    rich_help_panel="Commands" # Optional: for better help output with Rich
)

app.add_typer(mixnorm_ksigma.app, name="mixnorm_ksigma", help="Approximate RLCT of a pymc model")


if __name__ == "__main__":
    app()
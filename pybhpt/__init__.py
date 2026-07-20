from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("pybhpt")
except PackageNotFoundError:  # running from a source tree without an install
    __version__ = "unknown"

__all__ = ["__version__"]

"""Filesystem layout for the analysis scripts.

Raw downloads and processed intermediates are large and never stored in the repository. Point
``NK_DATA_DIR`` at a directory outside the repository (default: ``./data``) and the scripts use
``<NK_DATA_DIR>/raw`` and ``<NK_DATA_DIR>/processed``.
"""
import os

DATA_DIR = os.path.abspath(os.environ.get("NK_DATA_DIR", "data"))
RAW_DIR = os.path.join(DATA_DIR, "raw")
PROCESSED_DIR = os.path.join(DATA_DIR, "processed")

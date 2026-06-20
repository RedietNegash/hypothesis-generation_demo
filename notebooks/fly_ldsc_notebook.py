# -*- coding: utf-8 -*-
import marimo

__generated_with = "0.9.14"
app = marimo.App(width="medium")


@app.cell
def __():
    import marimo as mo
    import os
    import re
    import subprocess
    import pandas as pd
    import numpy as np
    from pathlib import Path
    import glob
    import concurrent.futures
    import multiprocessing
    return mo, os, re, subprocess, pd, np, Path, glob, concurrent, multiprocessing


if __name__ == "__main__":
    app.run()

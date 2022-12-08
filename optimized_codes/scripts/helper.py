from __future__ import annotations
import os
import numpy as np
import pickle
import click
import pathlib
import scanpy as sc

def write_anndata(anndata: sc.AnnData, dir: pathlib.Path, filename: str) -> None:
    """
    Writes scanpy AnnData object to disk in specified directory with filename
    """
    fullpath = pathlib.Path(dir, filename)
    click.echo(f"Writing AnnData object to {fullpath}")
    anndata.write(fullpath)

def create_color_palette(output_path: str | pathlib.Path, n_colors:int=10, overwrite:bool=False) -> None:
    """
    Create a custom color palette for plotting

    Args:
        output_path (str | pathlib.Path): _description_
        n_colors (int, optional): _description_. Defaults to 10.
        overwrite (bool, optional): _description_. Defaults to True.
    """
    my_colors = [np.random.rand(3,) for _ in range(n_colors)]
    
    if os.path.exists(output_path):
        if overwrite:
            print("Warning: overwriting existing color palette")
            os.remove(output_path)
        else:
            print("Warning: color palette already exists")
            print("If you want to overwrite existing palette, set overwrite=True")
            return
    
    with open(output_path, 'wb') as f:
       pickle.dump(my_colors, f)
       print("Color palette saved to:", os.path.abspath(output_path))
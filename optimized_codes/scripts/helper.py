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
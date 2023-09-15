import setuptools

with open("README.md", "r", encoding="utf-8") as fhand:
    long_description = fhand.read()

setuptools.setup(
    name="pcrglobwb_eval",
    version="0.0.1",
    author="Barry van Jaarsveld",
    author_email="a.s.vanjaarsveld@uu.nl",
    description=("A package to evaluate PCR-GLOBWB simulations."),
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=["xarray", "zarr", "geopandas", "hydromt", "pandas", "tqdm"],
    packages=setuptools.find_packages(),
    entry_points={
        "console_scripts": [
            "zview = zview.cli:main",
        ]
    }
)

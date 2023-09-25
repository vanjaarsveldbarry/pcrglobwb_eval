import setuptools


setuptools.setup(
    name="pcrglobwb_eval",
    version="0.0.1",
    author="Barry van Jaarsveld",
    author_email="a.s.vanjaarsveld@uu.nl",
    description=("A package to evaluate PCR-GLOBWB simulations."),
    long_description='test',
    long_description_content_type="text/markdown",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    packages=setuptools.find_packages(),
)
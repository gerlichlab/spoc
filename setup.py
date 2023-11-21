#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages


requirements = [
    "Click>=7.0",
    "pandas==1.3.0",
    "numpy>=1.19.0",
    "pandera>=0.9.0",
    "pydantic==1.9.1",
    "pyarrow==8.0.0",
    "hypothesis==6.52.4",
    "dask>=2023.3.0",
    "ipykernel",
    "matplotlib>=3.5.3",
    "seaborn>=0.11.2",
    "duckdb==0.9.1",
    "sparse==0.13.0",
    "numba>=0.57.0",
]

test_requirements = [
    "pytest>=3",
]

setup(
    author="Michael Mitter",
    author_email="michael.mitter@imba.oeaw.ac.at",
    python_requires=">=3.8",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
    description="Utility package for processing sister-sensitive PoreC data",
    entry_points={
        "console_scripts": [
            "spoc=spoc.cli:main",
        ],
    },
    install_requires=requirements,
    license="MIT license",
    include_package_data=True,
    keywords="spoc",
    name="spoc",
    packages=find_packages(include=["spoc", "spoc.*"]),
    test_suite="tests",
    tests_require=test_requirements,
    url="https://github.com/mittmich/spoc",
    version="0.1.1",
    zip_safe=False,
)

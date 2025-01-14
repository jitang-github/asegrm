from setuptools import setup, Extension
import numpy as np

# def readme():
#     with open("README.md") as f:
#         return f.read()

setup(
    name="asegrm",
    version="0.1",
    author="Ji Tang",
    author_email="jitang@usc.edu",
    description="Ancestry-specific Genetic Relationship Matrix",
    # long_description=readme(),
    long_description_content_type="text/markdown",
    # url="xxxx",
    # packages=["asegrm"],
    install_requires=[
        "pandas==1.5.3",
        "tskit==0.5.6",
        "tqdm==4.67.1",
    ],
    scripts=[
        "src/asegrm",
    ],
    ext_modules=[
        Extension("asegrm_matrix", ["src/asegrm_matrix.c"], include_dirs=[np.get_include()]),
    ],
    classifiers=[
        "Programming Language :: Python :: 3.9 :: CPython",
        "License :: OSI Approved :: MIT License",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    keywords="genetic-relationship-matrix admixture ancestry genealogy",
)

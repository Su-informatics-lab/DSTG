import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="DSTG",
    version="0.0.1",
    author="QSong",
    author_email="wasqqdyx@gmail.com",
    description=
    "Deconvoluting Spatial Transcriptomics data through Graph-based convolutional networks (DSTG)",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Su-informatics-lab/DSTG",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
)

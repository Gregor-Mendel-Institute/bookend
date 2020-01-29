import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="bookend-maschon0", # Replace with your own username
    version="0.0.1",
    author="Michael A. Schon",
    author_email="michael.schon@gmi.oeaw.ac.at",
    description="End-guided transcript assembler for short and long RNA-seq reads.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Gregor-Mendel-Institute/bookend",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    keywords='transcriptome assembler bioinformatics rna sequencing',
    install_requires=['cython', 'numpy'],
    python_requires='>=3.6',
)


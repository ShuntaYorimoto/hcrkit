from setuptools import setup, find_packages

setup(
    name="hcrkit",
    version="2.0.0",
    description="Automated pipeline for HCR (Hybridization Chain Reaction) probe design",
    author="Shunta Yorimoto",
    url="https://github.com/ShuntaYorimoto/hcrkit",
    packages=find_packages(),
    py_modules=['core'],
    scripts=['hcrkit.py'],
    install_requires=[
        'biopython>=1.78',
        'pandas>=1.3.0',
    ],
    python_requires='>=3.7',
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
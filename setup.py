from setuptools import setup, find_packages

setup(
    name="hcrkit",
    version="1.0",
    description="Automated pipeline for HCR (Hybridization Chain Reaction) probe design with customizable parameters",
    author="Shunta Yorimoto",
    url="https://github.com/ShuntaYorimoto/hcrkit",
    packages=find_packages(),
    entry_points={
        'console_scripts':[
            'hcrkit=hcrkit:main',
            'extract_target_ids=extract_target_ids:main'
        ]
    },
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
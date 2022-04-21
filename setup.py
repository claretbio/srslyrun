import setuptools
import subprocess
import os

with open("README.md", "r") as fh:
    long_description = fh.read()


subprocess.run(['tar', '-czvf', 'srsly/workflow.tar.gz', '-C', '/mnt/svetlana/sandbox/camille/DL-302/srsly-software/', 'workflow'])


setuptools.setup(
    name="srslyrun",
    version="0.1.0",
    author="Camille Schwartz",
    author_email="camille@claretbio.com",
    description="Not your grandma's NGS analysis - software for analyzing FASTQs from SRSLY libraries",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/claretbio/srslyrun",
    packages=setuptools.find_packages(),
    test_suite='nose.collector',
    tests_require=['nose', 'coverage'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
            'snakemake==7.3.1',
            'umi_tools',
            'pyyaml',
            'srslyumi',
            'pysam'
        ],
    package_data={'srsly': ['workflow.tar.gz']},
    include_package_data=True,
    entry_points={"console_scripts": ["srsly=srsly.cli:main"]},
    python_requires=">=3.7",
)

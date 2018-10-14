#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()

REQUIREMENTS = [
    'mavis>=2.1.4'
]

SETUP_REQUIREMENTS = ['pytest-runner', ]

TEST_REQUIREMENTS = [
    'pip==9.0.1',
    'bumpversion==0.5.3',
    'wheel==0.30.0',
    'watchdog==0.8.3',
    'flake8==3.5.0',
    'tox==2.9.1',
    'coverage==4.5.1',
    'Sphinx==1.7.1',
    'twine==1.10.0',
    'pytest>=3.4.2',
    'pytest-runner>=2.11.1',
    'pytest-cov>=2.6.0'
]

setup(
    author="Morgan Bye",
    author_email='morgan@morganbye.com',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    description='Simulating transcriptome structural variants',
    entry_points={
        'console_scripts': [
            'rna_sv_simulator=rna_sv_simulator.cli:main',
        ],
    },
    extras_require={
        'dev': TEST_REQUIREMENTS 
    },
    install_requires=REQUIREMENTS,
    license="GNU General Public License v3",
    long_description=readme,
    include_package_data=True,
    keywords='rna_sv_simulator',
    name='rna_sv_simulator',
    packages=find_packages(include=['rna_sv_simulator']),
    setup_requires=SETUP_REQUIREMENTS,
    test_suite='tests',
    tests_require=TEST_REQUIREMENTS,
    url='https://github.com/morganbye/rna_sv_simulator',
    version='0.1.0',
    zip_safe=False,
)

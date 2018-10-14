#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()

requirements = ['yaml']

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest', ]

setup(
    author="Cara Reisle",
    author_email='creisle@bcgsc.ca',
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
        'dev': [
            'pip==9.0.1',
            'bumpversion==0.5.3',
            'wheel==0.30.0',
            'watchdog==0.8.3',
            'flake8==3.5.0',
            'tox==2.9.1',
            'coverage==4.5.1',
            'Sphinx==1.7.1',
            'twine==1.10.0',
            'pytest==3.4.2',
            'pytest-runner==2.11.1',
            'pytest-fs',
            'hypothesis',
        ]
    },
    install_requires=requirements,
    license="GNU General Public License v3",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='rna_sv_simulator',
    name='rna_sv_simulator',
    packages=find_packages(include=['rna_sv_simulator']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/creisle/rna_sv_simulator',
    version='0.1.0',
    zip_safe=False,
)

import os
import unittest

from setuptools import find_packages, setup
from setuptools.command.build_py import build_py

long_description = ""
if os.path.exists("README.md"):
    with open("README.md", "rt") as fh:
        long_description = fh.read()

with open("requirements.txt", "rt") as f:
    requirements = [r.strip() for r in f.readlines()]

class PublishCommand(build_py):
    """Publish package to PyPI"""
    def run(self):
        os.system("rm -rf dist")
        os.system("python3 setup.py sdist"
                  "&& python3 setup.py bdist_wheel"
                  "&& python3 -m twine upload dist/*whl dist/*gz")


def test_suite():
    """Discover unittests"""
    test_loader = unittest.TestLoader()
    test_suite = test_loader.discover('annotation_utils', pattern='*tests.py')
    return test_suite


setup(
    name='annotation-utils',
    version="0.1.4",
    description="Misc. utilities for downloading and working with various gene annotations, the HPO ontology, etc",
    install_requires=requirements,
    cmdclass={
        'publish': PublishCommand,
    },
    entry_points = {
        'console_scripts': [
            'hpo_lookup = annotation_utils.hpo_lookup:main',
        ],
    },
    long_description_content_type="text/markdown",
    long_description=long_description,
    packages=["annotation_utils"],
    include_package_data=True,
    python_requires=">=3.7",
    license="MIT",
    keywords='',
    test_suite="setup.test_suite",
    #tests_require=["mock"],
    url='https://github.com/bw2/annotation-utils',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: Implementation :: CPython',
        'Programming Language :: Python :: Implementation :: PyPy',
    ],
)

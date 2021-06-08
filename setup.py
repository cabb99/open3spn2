"""
Open3SPN2
An Implementation of the 3SPN.2 and 3SPN.2C coarse-grained molecular models of DNA in OpenMM.
"""
import sys
from setuptools import setup, find_packages
import versioneer

print(f"Version {versioneer.get_version()}")
short_description = __doc__.split("\n")

try:
    with open("README.md", "r") as handle:
        long_description = handle.read()
except:
    long_description = "\n".join(short_description[2:])

setup(
      name = 'open3SPN2',
      packages = ['open3SPN2'],
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
      license='MIT',
      description = short_description[0],
      long_description=long_description,
      long_description_content_type="text/markdown",
      author = 'Carlos Bueno', 
      author_email = 'carlos.bueno@rice.edu', 
      url = 'https://github.com/cabb99/open3spn2', 
      download_url = 'https://github.com/cabb99/open3spn2/archive/0.2.0.tar.gz',
      keywords = ['dna', 'forcefield', 'openmm'],
      package_data={'open3SPN2': ['3SPN2.conf','3SPN2.xml']},
      include_package_data=True,
      install_requires=[
              'biopython',
              'pandas',
              'numpy',
              'scipy',
              'mdtraj',
              'openmm',
              'pdbfixer',
              'nose',
          ],
      classifiers=[
        'Development Status :: 3 - Alpha', 
        'Intended Audience :: Science/Research', 
        'Topic :: Scientific/Engineering :: Physics', 
        'License :: OSI Approved :: MIT License', 
        'Programming Language :: Python :: 3', 
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
      ],
    )

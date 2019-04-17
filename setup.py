from setuptools import setup, find_packages


with open('README.md') as f:
  readme = f.read()

with open('LICENSE') as f:
  license = f.read()

setup(
  name='EverFORC',
  version='0.0.4',
  description='Script to generate sequence files for an Evercool II to perform FORC measurements, as well as to analyze them',
  url='https://github.com/DomiDre/EverFORC',
  author='Dominique Dresen',
  author_email='dominique.dresen@uni-koeln.de',
  license=license,
  long_description=readme,
  install_requires=[
    'numpy',
    'matplotlib'
  ],
  python_requires='>2.7, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*',
  platforms=['Linux'],
  package_dir={'EverFORC': 'EverFORC'},
  packages=find_packages(
    exclude=(
      '_build',
      'docs',
      '_static',
      '_templates'
      'tests',
      'examples'
      )
  ),
  keywords='ppms evercool FORC magnetism'
)
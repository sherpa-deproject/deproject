from setuptools import setup, find_packages

with open('README.rst', 'rt') as fh:
    long_description = fh.read()

setup(name='deproject',
      version='0.1.3',
      license='BSD',
      description='Sherpa deprojection package (X-ray analysis of Galaxy Clusters, Groups, and Galaxies)',
      keywords='deprojection xray 3d 2d plasma Astrophysics onion',
      long_description=long_description,
      long_description_content_type='text/x-rst',
      author='Tom Aldcroft',
      author_email='aldcroft@cfa.harvard.edu',
      # To be changed to read-the-docs
      url='http://cxc.harvard.edu/contrib/deproject/',

      packages=find_packages('src'),
      package_dir={'': 'src'},

      # QUS: Should the test data files be included?

      # NOTE:
      #  CIAO 4.11 comes with NumPy 1.12.1, but AstroPy 3.1 requires
      #  NumPy >= 1.13.1, so add the restriction here.
      #
      #  This is restrictive for standalone Sherpa users, but it
      #  requires a *lot* of effort to build the needed XSPEC model
      #  support into Sherpa, so assume that the majority of use
      #  cases will be installing into CIAO.
      #
      # SciPy is needed by astropy.cosmology
      #
      setup_requires=['pytest-runner'],
      tests_require=['pytest'],
      install_requires=['sherpa', 'astropy<3.1', 'scipy'],

      python_requires='~=3.5',

      classifiers=[
          'Development Status :: 4 - Beta',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: BSD License',
          'Programming Language :: Python :: 3.5',
          'Programming Language :: Python :: 3.6',
          'Programming Language :: Python :: 3.7',
          'Topic :: Scientific/Engineering :: Astronomy',
          'Topic :: Scientific/Engineering :: Physics'
          ],

      # I think this should be zip_safe, but have not tried it
      # zip_safe=False
      )

from setuptools import setup, find_packages

with open('README.rst', 'rt') as fh:
    long_description = fh.read()

setup(name='deproject',
      version='0.2.1',
      license='BSD',
      description='Sherpa deprojection package (X-ray analysis of Galaxy Clusters, Groups, and Galaxies)',
      keywords='deprojection xray 3d 2d plasma Astrophysics onion',
      long_description=long_description,
      long_description_content_type='text/x-rst',
      author='Douglas Burke, Tom Aldcroft',
      author_email='dburke.gw@gmail.com',

      url='https://deproject.readthedocs.io/',
      project_urls={
          'Documentation': 'https://deproject.readthedocs.io/',
          'Source Code': 'https://github.com/sherpa-deproject/deproject/',
          'Issues': 'https://github.com/sherpa-deproject/deproject/issues/',
          },

      packages=find_packages('src'),
      package_dir={'': 'src'},

      setup_requires=['pytest-runner'],
      tests_require=['pytest'],
      install_requires=['sherpa', 'astropy', 'scipy'],

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

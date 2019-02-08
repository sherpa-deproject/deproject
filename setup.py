from setuptools import setup

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
      url='http://cxc.harvard.edu/contrib/deproject/',
      py_modules=['deproject', 'specstack', 'cosmocalc'],
      install_requires=['sherpa'],
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

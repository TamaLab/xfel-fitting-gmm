from setuptools import setup

setup(name='AFM',
      version='0.1',
      description='Sampling against AFM/XFEL data',
      long_description='Sampling against AFM/XFEL data',
      url='',
      license='Apache2.0',
      packages=['AFM'], 
      include_package_data=True, 
      test_suite='nose.collector',
      tests_require=['nose'],
      package_data={'AFM': ['clib/*.so']},
      scripts=['bin/afmJoke','bin/afmEmulator','bin/afmSampler',
               'bin/volCompare12','bin/imgCompare12','bin/overlap12',
               'bin/mapCompare12','bin/makeMap','bin/fitGmmRigid12',
               'bin/trjReader'], 
      entry_points = {
        'console_scripts': ['afmJoke=AFM.command_line:main'],
      },
      zip_safe=False)

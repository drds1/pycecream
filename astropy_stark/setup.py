from setuptools import setup

#pip install .
#python setup.py sdist
#twine upload dist/astropy_stark-1.0.1.tar.gz
setup(name='astropy_stark',
      version='1.0.1',
      description='random-custom astronomy functions including making fake accretion disk light curves',
      url='https://github.com/dstarkey23/academic_projects_public',
      author='dstarkey23',
      author_email='ds207@st-andrews.ac.uk',
      license='MIT',
      packages=['astropy_stark'],
      #packages=['fish_forecast'],
      install_requires=[
      'scikit-learn',
      'scipy',
      'statsmodels',
      'pandas',
      'numpy',
      'matplotlib',
      'glob3'
      #'prediction_functions',
      ],
      zip_safe=False)
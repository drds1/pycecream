from setuptools import setup

#pip install .
#python setup.py sdist
#twine upload dist/astropy_stark-1.1.9.tar.gz
setup(name='astropy_stark',
      version='1.1.9',
      description='random-custom astronomy functions including making fake accretion disk light curves',
      url='https://github.com/dstarkey23/pycecream',
      author='dstarkey23',
      author_email='ds207@st-andrews.ac.uk',
      license='MIT',
      packages=['astropy_stark'],
      install_requires=[
      'scikit-learn',
      'scipy',
      'statsmodels',
      'pandas',
      'numpy',
      'matplotlib',
      'glob3',
      ],
      zip_safe=False)

#jupyter nbconvert --to Markdown ./pycecream_test/test_pycecream.ipynb
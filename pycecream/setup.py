from setuptools import setup

#uload to pip
#pip install .
#python setup.py sdist
#twine upload dist/pycecream-1.0.0.tar.gz

setup(name='pycecream',
      version='1.0.1',
      description='python implementation of the cream accretion disc fitting code '
                  'https://academic.oup.com/mnras/article-abstract/456/2/1960/1066664?redirectedFrom=PDF',
      url='https://github.com/dstarkey23/academic_projects_public',
      author='dstarkey23',
      author_email='ds207@st-andrews.ac.uk',
      license='MIT',
      packages=['pycecream'],
      #packages=['fish_forecast'],
      install_requires=[
      #'scikit-learn',
      'pandas',
      'numpy',
      'matplotlib',
      'scipy',
      'astropy_stark',
      'glob3'
      #statsmodels',
      #'pandas'
      #'prediction_functions',
      ],
      zip_safe=False)
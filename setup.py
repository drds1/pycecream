from setuptools import setup

#upload to pip
#pip install .
#python setup.py sdist
#twine upload dist/pycecream-1.5.5.tar.gz

setup(name='pycecream',
      version='1.5.5',
      description='python implementation of the cream accretion disc fitting code '
                  'https://academic.oup.com/mnras/article-abstract/456/2/1960/1066664?redirectedFrom=PDF'
      ,
      long_description= 'add pycecream.dream light curve merging feature'
      ,url='https://github.com/dstarkey23/pycecream',
      author='dstarkey23',
      author_email='ds207@st-andrews.ac.uk',
      license='MIT',
      packages=['pycecream'],
      package_data={'': ['creaminpar.par','cream_f90.f90']},
      install_requires=[
      'pandas==2.2.0',
      'numpy==1.26.4',
      'matplotlib==3.8.2',
      'scipy==1.12.0',
      'corner==2.2.2',
      'seaborn==0.13.2'
      ],
extras_require={
        'tests': [
            'nose2==0.9.1',
            'pre-commit==1.20.0',
            'flake8==3.7.9',
            'pydoc-markdown==2.0.4',
            'tabulate==0.8.5',
            'six==1.12.0'
        ]
    },
      zip_safe=False)

#jupyter nbconvert --to Markdown ./pycecream_test/test_pycecream.ipynb
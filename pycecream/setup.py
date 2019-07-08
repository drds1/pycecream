from setuptools import setup

#upload to pip
#pip install .
#python setup.py sdist
#twine upload dist/pycecream-1.2.2.tar.gz

setup(name='pycecream',
      version='1.2.2',
      description='python implementation of the cream accretion disc fitting code '
                  'https://academic.oup.com/mnras/article-abstract/456/2/1960/1066664?redirectedFrom=PDF'
      ,
      long_description= 'bug fix to return the correct rescaled error bars in output_merged_data dictionary in the get_light_curve_fits function'
      'and scaling parameters',
      url='https://github.com/dstarkey23/pycecream',
      author='dstarkey23',
      author_email='ds207@st-andrews.ac.uk',
      license='MIT',
      packages=['pycecream'],
      package_data={'': ['creaminpar.par','cream_f90.f90']},
      install_requires=[
      'pandas',
      'numpy',
      'matplotlib',
      'scipy',
      'astropy_stark',
      'glob3'
      ],
      zip_safe=False)

#jupyter nbconvert --to Markdown ./pycecream_test/test_pycecream.ipynb
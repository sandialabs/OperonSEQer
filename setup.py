import os
from setuptools import setup

install_requires = ['pandas==1.3','numpy==1.22','configargparse','seaborn==0.11.2','scipy==1.7.1','scikit-learn==1.0','matplotlib==3.4.3', 'xgboost==1.4.2', 'pyinstaller==3.6', 'libomp==11.1.0']
#install_requires=[]

setup(
    name="OperonSEQer", 
    version='1.0', 
    author = 'Raga Krishnakumar, Anne Ruffing',
    description = 'RNA-sequencing data and multiple ML models used to determine operon pairs in prokaryotes',
    install_requires=install_requires,
    license='MPL 2.0'
)

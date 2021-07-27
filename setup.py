import os
from setuptools import setup

install_requires = ['pandas==0.23.4','numpy==1.16','configargparse','seaborn==0.9.0','scipy==1.1.0','scikit-learn==0.20.1','matplotlib', 'xgboost==1.2.0']
#install_requires=[]

setup(
    name="OperonSEQer", 
    version='1.0', 
    author = 'Raga Krishnakumar, Anne Ruffing',
    description = 'RNA-sequencing data and multiple ML models used to determine operon pairs in prokaryotes',
    install_requires=install_requires,
    license='MPL 2.0'
)
from setuptools import setup, find_packages

setup(
    name='map3C',
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    #setup_requires=['setuptools'],
    author='Joseph Galasso',
    author_email='jgalasso1@g.ucla.edu',
    packages=find_packages(),
    description='A package for mapping 3C/Hi-C data',
    long_description=open('README.md').read(),
    license='MIT',
    long_description_content_type='text/markdown',
    url='https://github.com/solveforj/map3C',
    include_package_data=True,
    #install_requires=['numpy', 'scipy', 'scikit-learn', 'h5py', 'joblib', 'clodius',
    #                  'tables', 'cooler', 'pandas', 'statsmodels', 'rpy2', 'anndata', 'xarray', 'zarr', 'numcodecs'],
    package_data={
        '': ['*.txt', '*.Snakefile', '*.smk']
    },
    entry_points={
        'console_scripts': ['map3C=map3C.__main__:main'],
    }
)
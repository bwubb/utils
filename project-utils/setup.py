from setuptools import setup,find_packages

setup(
    name='project-utils',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        'pyyaml',
    ],
    entry_points={
        'console_scripts':[
            'project-utils=project_utils.cli:main',
        ],
    },
) 
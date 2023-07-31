from setuptools import setup

setup(
    name='GoCa evaluation',
    version='1.0.0',
    author='Luis Walter',
    author_email='luis.walter@tum.de',
    description='Evaluate structure-based coarse-grained protein model simulations',
    scripts=['analysis.py', 'gromacs_analysis.py'],
    install_requires=['numpy', 'scipy', 'mdtraj', 'networkx', 'matplotlib', 'pandas'],
)
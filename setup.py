from setuptools import setup, find_packages
setup(
    name="TMB",
    version="0.1",
    packages=find_packages(exclude=["tests"]),
    scripts=['*.py', 'main.py'],
)
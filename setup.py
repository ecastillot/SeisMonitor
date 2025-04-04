from setuptools import setup, find_packages
import pathlib

import pkg_resources
import setuptools
import codecs
import os
import SeisMonitor
# here = os.path.abspath(os.path.dirname(__file__))

# with codecs.open(os.path.join(here, "README.md"), encoding="utf-8") as fh:
#     long_description = "\n" + fh.read()

VERSION = SeisMonitor.__version__
DESCRIPTION = 'To monitor seismic activity'
LONG_DESCRIPTION = 'A package that allows to monitor the seismic activity through main steps in the monitoring workflow: earthquake detection and phase picking -> phase associator -> earthquake locator -> magnitude estimation.'

req_path = os.path.join(os.path.dirname(__file__),"requirements.txt")
with pathlib.Path('requirements.txt').open() as requirements_txt:
    install_requires = [
        str(requirement)
        for requirement
        in pkg_resources.parse_requirements(requirements_txt)
    ]

# Setting up
setup(
    name="seismonitor",
    version=VERSION,
    author="ecastillot (Emmanuel Castillo)",
    author_email="<castillo.280997@gmail.com>",
    url="https://github.com/ecastillot/SeisMonitor",
    description=DESCRIPTION,
    long_description_content_type="text/markdown",
    long_description=LONG_DESCRIPTION,
    packages=find_packages(),
    install_requires=install_requires,
    keywords=['python', "seismonitor","earthquakes","seismology"],
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Developers",
        "Programming Language :: Python :: 3",
        "Operating System :: Unix",
    ],
    python_requires='>=3.8'
)

# python setup.py sdist bdist_wheel
# twine upload dist/*
# python -m twine upload -u __token__ -p [unique_token] dist/*

#pip install twine
# pip install sphinx
#pip install sphinx-rtd-theme
# https://www.youtube.com/watch?v=BWIrhgCAae0&ab_channel=MariattaWijaya
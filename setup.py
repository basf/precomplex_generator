"""
# Installer for the precomplex generator scripts
"""

import os

import setuptools
# numpy is needed for 'f2py' code translation/compilation
from numpy.distutils.core import Extension, setup


def main():
    """Main function"""

    # packages
    packages = setuptools.find_packages(exclude=["*.tests"])

    # Entry points
    # A script will be automatically created for each entry point that calls the function from the given module
    entry_points = {
        "console_scripts": [
            "prec_gen = precomplex_generator.main:main",
        ],
    }

    #
    # Package dependencies
    #
    # dependencies = open(
    #     os.path.join(os.path.dirname(__file__), "requirements.txt")
    # ).readlines()

    # dependencies_setup = [
    #     "setuptools-git >= 0.4",
    # ]

    #
    # setup
    #
    setup(
        name="precomplex_generator",
        description="Generator of a precomplex structure, i.e. a suited input structure for single-ended reaction path optimization.",
        author="Julian Geiger, Maike Bergeler",
        author_email="jgeiger@iciq.es, maike.bergeler@basf.com",
        license="BASF proprietary",
        long_description="",
        url="https://github.com/basf/precomplex_generator.git",
        packages=packages,
        entry_points=entry_points,
        # ext_modules=extmodules,
        # install_requires=dependencies,
        # setup_requires=dependencies_setup,
        zip_safe=False,
    )


if __name__ == "__main__":
    main()

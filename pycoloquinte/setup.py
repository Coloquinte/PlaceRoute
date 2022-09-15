from setuptools import setup

with open('../README.md') as f:
    long_description = f.read()

setup(
    name="coloquinte",
    version="0.0.1",
    author="Gabriel Gouvine",
    author_email="gabriel.gouvine_git@m4x.org",
    url="https://github.com/Coloquinte/PlaceRoute",
    description="Python interface for the Coloquinte VLSI placer",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="MIT",
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Developers",
        "Topic :: Scientific/Engineering :: Electronic Design Automation (EDA)",
    ],
    py_modules=["coloquinte"],
    dependencies=["coloquinte_pybind"],
    include_package_data=True,
    python_requires=">=3.6",
)

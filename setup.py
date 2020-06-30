import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="SS_MTI",  # Replace with your own username
    version="0.0.1",
    author="Nienke Brinkman",
    author_email="nienke.brinkman@erdw.ethz.ch",
    description="A single-station source inversion package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/nienkebrinkman/SS_MTI/",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
)

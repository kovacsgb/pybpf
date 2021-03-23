import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pybpf", # Replace with your own username
    version="0.9.0",
    author="Gábor Kovács",
    author_email="kovacs.gabor@csfk.org",
    description="Interface package for working with Budapest-Florida-code outputs",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/kovacsgb/pybpf",
    project_urls={
        "Bug Tracker": "https://github.com/kovacsgb/pybpf/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    packages=['pybpf'],
    python_requires=">=3.4",
    install_requires = ['numpy','astropy']
)
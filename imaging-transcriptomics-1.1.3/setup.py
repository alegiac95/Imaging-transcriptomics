from setuptools import setup, find_packages
from pathlib import Path
from glob import glob


def read_long_description():
    """Read the README file and use it as long description in Pypi.

    :return: contents of the readme file.
    """
    with open("README.md", "r") as f:
        return f.read()


def get_version_number():
    """Get the version of the package from the __init__ file.

    :return version: version of the package.
    """
    with open("imaging_transcriptomics/__init__.py", "r") as f:
        for line in f:
            if '__version__' in line:
                _, version, _ = line.split('"')
                break
        return version


def get_requirements():
    """Get the requirements for the installation."""
    requirements = []
    with open('requirements.txt', 'r') as f:
        for line in f:
            if line[0] == '#':
                continue
            requirements.append(line.rstrip('\n'))
    return requirements


setup(name="imaging-transcriptomics",
      author="Alessio Giacomel, Daniel Martins",
      author_email="alessio.giacomel@kcl.ac.uk , daniel.martins@kcl.ac.uk",
      version=get_version_number(),
      url="https://github.com/alegiac95/Imaging-transcriptomics",
      description="A package to perform imaging transcriptomics analysis on a "
                  "neuroimaging brain scan.",
      license="GPLv3",
      long_description=read_long_description(),
      long_description_content_type="text/markdown",
      classifiers=["Intended Audience :: Healthcare Industry",
                   "Intended Audience :: Science/Research",
                   "Topic :: Scientific/Engineering :: Image Processing",
                   "Topic :: Scientific/Engineering :: Medical Science Apps.",
                   "Development Status :: 4 - Beta",
                   "Programming Language :: Python :: 3",
                   "Programming Language :: Python :: 3.6",
                   ],
      keywords="Image analysis, Neuroimaging, Imaging Transcriptomics, "
               "Medical Imaging, Research, Multimodal Imaging",
      install_requires=get_requirements(),
      packages=["imaging_transcriptomics"],
      include_package_data=True,
      package_data={"imaging_transcriptomics": ["*.yaml", "*.nii.gz",
                                                "*.csv", "*.txt", "*.annot"]},
      entry_points={
          "console_scripts": [
              "imagingtranscriptomics=imaging_transcriptomics.script"
              ".imagingtranscriptomics:main",
              "imt_gsea=imaging_transcriptomics.script.imt_gsea:main",
          ]
      },
      project_urls={
        "Source": "https://github.com/alegiac95/Imaging-transcriptomics",
        "Bug Reports": "https://github.com/alegiac95/Imaging-transcriptomics/issues",
        "Documentation": "https://imaging-transcriptomics.rtfd.io/"
      },
      )

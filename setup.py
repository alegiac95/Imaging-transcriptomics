from setuptools import setup, find_packages


def read_long_description():
    """Read the README file and use it as long description in Pypi.

    :return: contents of the readme file.
    """
    with open("README.rst", "r") as f:
        return f.read()


def get_version_number():
    with open("imaging_transcriptomics/__init__.py", "r") as f:
        line = f.readline()
        version = line.split(" = ")[1]
        return version


setup(name="imaging-transcriptomics",
      author="Daniel Martins, MD, PhD; Alessio Giacomel",
      author_email=["daniel.martins@kcl.ac.uk", "alessio.giacomel@kcl.ac.uk"],
      version=get_version_number(),
      description="Analyse an imaging scan with imaging transcriptomics",
      long_description=read_long_description(),
      classifiers=["Intended Audience :: Healthcare Industry",
                   "Intended Audience :: Science/Research",
                   "Topic :: Scientific/Engineering :: Image Processing",
                   "Topic :: Scientific/Engineering :: Medical Sciences App",
                   "Development Status :: 4 - Beta",
                   "Programming Language :: Python :: 3",
                   "Programming Language :: Python :: 3.6",
                   ],
      keywords="Image analysis, Neuroimaging, Imaging Transcriptomics, Medical Imaging, Research, Multimodal Imaging",
      packages=find_packages(),
      entry_points={"console_scripts": ["imaging-transcriptomics = scripts.imagingtranscriptomics:main"]}
      )


import argparse

import imaging_transcriptomics


def get_args():
    """Parse the inputs to run the virtual histology script.
    """
    DESCRIPTION = """Perform imaging transcriptomics analysis on a neuroimaging scan. """
    EPILOG = """Check your results in the specified folder or in the file path of the input scan, if you have not 
    specified an output path. If you used this software in your research please cite: 
    * Imaging transcriptomics: **Convergent cellular, transcriptomic, and molecular neuroimaging signatures in the 
    healthy adult human brain.**
    *Daniel Martins, Alessio Giacomel, Steven CR Williams, Federico E Turkheimer, Ottavia Dipasquale, Mattia Veronese, 
    PET templates working group*; bioRxiv 2021.06.18.448872; doi: `https://doi.org/10.1101/2021.06.18.448872 
    <https://doi.org/10.1101/2021.06.18.448872>`_ """

    parser = argparse.ArgumentParser(description=DESCRIPTION,
                                     epilog=EPILOG)

    parser.add_argument("-i", "--input",
                        type=str,
                        help=("Input imaging file in NIfTI format (.nii, .nii.gz).\n"
                              "The input file is expected to have the same matrix size as the atlas used (182x218x182),"
                              "if the input image has different matrix size this can be resliced to match the"
                              "resolution of the MNI152 1mm template matrix size (e.g. the one provided with FSL)."),
                        required=True)
    parser.add_argument("-o", "--out",
                        type=str,
                        help="Path where to save the output, if not specified the path of the path of the input scan "
                             "is used.",
                        required=False)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-n", "--ncomp",
                       type=int,
                       help="Number of PLS components to use. The number of components has to be between 1 and 15.")
    group.add_argument("-v", "--variance",
                       type=int,
                       help="""Variance explained by the components. The variance input should be between 10 and 
                       100, and the program will select the number of components that explain a variance closest to 
                       the desired (with the explained variance used as a minimum). """)
    return parser.parse_args()


def main():
    inputs = get_args()
    # TODO: use multiprocessors if the input is a list of strings.

    data_to_analyse = imaging_transcriptomics.inputs.extract_average(
        imaging_transcriptomics.inputs.read_scan(inputs.input)
    )
    initial_dict = {
        "variance": inputs.variance,
        "n_components": inputs.ncomp
    }
    analysis = imaging_transcriptomics.ImagingTranscriptomics(data_to_analyse, **initial_dict)
    analysis.run()


if __name__ == '__main__':
    main()

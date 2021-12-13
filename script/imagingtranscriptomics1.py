#!/usr/bin/env python3

import argparse
from pathlib import Path
import imaging_transcriptomics as imt


def get_command_line_args():
    """
    Parse the command line and extract the relevant arguments to run the
    script.

    :return:
    """
    DESCRIPTION = ""
    EPILOG = ""

    parser = argparse.ArgumentParser(description=DESCRIPTION,
                                     epilog=EPILOG)
    # -- REQUIRED INPUTS --
    parser.add_argument("-i", "--input", type=str,
                        required=True,
                        help="")
    # Method to run the script
    method = parser.add_mutually_exclusive_group(required=True)
    method.add_argument("-c", "-corr", action="store_true", help="")
    method.add_argument("-p", "--pls", action="store_true", help="")


def main():
    """
    Main function to run the script.
    """
    parsed_inputs = get_command_line_args()
    # Directories and file checks
    in_file = parsed_inputs.input
    in_file = Path(in_file)
    if not in_file.exists():
        raise FileNotFoundError(f"{in_file} does not exist.")
    if in_file.suffix == ".txt":
        pass  # create the class from the @from_file method
    elif in_file.suffixes[-1] + in_file.suffixes[-2] == ".nii.gz" or \
            in_file.suffixes[-1] == ".nii":
        pass  # create the class from the @from_scan method


if __name__=="__main__":
    main()

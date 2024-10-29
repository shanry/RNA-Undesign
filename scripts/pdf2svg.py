#!/usr/bin/env python3
"""
This script converts PDF files to SVG format using Inkscape and optionally optimizes the SVG files using Scour.

Functions:
    process_file(args):
        Converts a single PDF file to SVG format using Inkscape.
        Args:
            args (tuple): A tuple containing the filename, input folder, and output folder.
    
    convert_pdf_to_svg(input_folder, output_folder):
        Converts all PDF files in the input folder to SVG format using multiprocessing.
        Args:
            input_folder (str): The folder containing the PDF files to be converted.
            output_folder (str): The folder where the converted SVG files will be saved.
    
    main():
        The main function that parses command-line arguments and initiates the PDF to SVG conversion process.
"""
import os
import subprocess
import sys
from multiprocessing import Pool, cpu_count

def process_file(args):
    filename, input_folder, output_folder = args
    if filename.endswith(".pdf"):
        input_pdf_path = os.path.join(input_folder, filename)
        output_svg_path = os.path.join(output_folder, filename.replace(".pdf", ".svg"))

        # Convert PDF to SVG using pdf2svg
        print(f"Converting {filename} to SVG...")
        # cmd_list = ["inkscape", "--without-gui", "--file", input_pdf_path, "--export-plain-svg", output_svg_path]
        cmd_list = [
            "inkscape",
            input_pdf_path,
            "--export-type=svg",
            "--export-filename", output_svg_path
        ]
        print(" ".join(cmd_list))
        subprocess.run(cmd_list)

        # # Use Scour to clean and crop the SVG
        # cropped_svg_path = os.path.join(output_folder, filename.replace(".pdf", "_cropped.svg"))
        # print(f"Cropping and optimizing {output_svg_path}...")
        # subprocess.run(["scour", "-i", output_svg_path, "-o", cropped_svg_path, "--enable-viewboxing"])

        # # Remove the uncropped SVG file to keep only the final cropped one
        # os.remove(output_svg_path)
        # os.rename(cropped_svg_path, output_svg_path)
        # print(f"Final cropped SVG saved as {output_svg_path}")

def convert_pdf_to_svg(input_folder, output_folder):
    # Create output folder if it doesn't exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Get a list of all PDF files in the input folder
    pdf_files = [(filename, input_folder, output_folder) for filename in os.listdir(input_folder) if filename.endswith(".pdf")]

    # Use multiprocessing to process files in parallel
    with Pool(cpu_count()) as pool:
        pool.map(process_file, pdf_files)

def main():
    if len(sys.argv) != 2:
        print("Usage: python convert_pdfs_to_svgs.py <input_folder>")
        sys.exit(1)

    input_folder = sys.argv[1]
    output_folder = os.path.join(input_folder, "svgs")

    # Convert all PDFs to SVG and crop them using Scour
    convert_pdf_to_svg(input_folder, output_folder)
    print(f"All PDF files converted and saved as SVG in {output_folder}.")

if __name__ == "__main__":
    main()
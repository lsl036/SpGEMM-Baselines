#!/usr/bin/env python3
import re
import csv
import sys
import argparse

def parse_log_file(input_file):
    """
    Parse the input log file to extract dataset names and times.
    Returns a list of (dataset, time) tuples.
    """
    results = []

    # Updated pattern to handle dataset names containing hyphens
    # This will match everything between "DATASET: " and "-RCM"
    dataset_pattern = re.compile(r'~+<DATASET: (.*?)>~+')
    time_pattern = re.compile(r'mtxfile_reorder: ([\d.]+) seconds')

    current_dataset = None

    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()

            # Check if the line contains dataset information
            dataset_match = dataset_pattern.search(line)
            if dataset_match:
                current_dataset = dataset_match.group(1)
                continue

            # Check if the line contains time information and we have a dataset
            time_match = time_pattern.search(line)
            if time_match and current_dataset:
                time_value = float(time_match.group(1))
                results.append((current_dataset, time_value))
                current_dataset = None  # Reset current dataset after recording a result

    return results

def save_to_csv(data, output_file):
    """
    Save the extracted data to a CSV file.
    """
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Dataset', 'Time (seconds)'])
        writer.writerows(data)

    print(f"Data successfully written to {output_file}")
    print(f"Extracted {len(data)} dataset-time pairs")

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Parse log file for dataset and time information')
    parser.add_argument('input_file', help='Path to the input log file')
    parser.add_argument('output_file', help='Path to the output CSV file')

    args = parser.parse_args()

    try:
        # Parse the log file
        results = parse_log_file(args.input_file)

        if not results:
            print("No dataset-time pairs found in the input file.")
            return

        # Save results to CSV
        save_to_csv(results, args.output_file)

    except FileNotFoundError:
        print(f"Error: File '{args.input_file}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
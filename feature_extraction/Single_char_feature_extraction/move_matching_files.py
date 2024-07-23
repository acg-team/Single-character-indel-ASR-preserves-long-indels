import os, sys
import shutil


if __name__ == '__main__':
    os.chdir(os.path.dirname(sys.argv[0]))
    # Step 1: Load the IDs from the first text file into a list
    ids_file = '../../mammals_data/00_agreeing_topology_IDs.txt'
    with open(ids_file, 'r') as file:
        ids = [line.strip() for line in file]

    # Step 2: Search the directory for files containing those IDs as part of their filenames
    source_directory = 'C:\\Users\\jowkar\\Documents\\Python scripts\\Single_char_indel_submision\\mammals_data\\output\\'
    destination_directory = 'C:\\Users\\jowkar\\Documents\\Python scripts\\Single_char_indel_submision\\mammals_data\\new_output\\'

    # Create the destination directory if it doesn't exist
    os.makedirs(destination_directory, exist_ok=True)

    selected_files = []

    for root, dirs, files in os.walk(source_directory):
        for file in files:
            for id in ids:
                if id in file:
                    full_file_path = os.path.join(root, file)
                    selected_files.append(full_file_path)
                    break  # Found a matching ID, no need to check other IDs for this file

    # Step 3: Move the selected files to the new directory
    for file in selected_files:
        shutil.move(file, destination_directory)
        print(f"Moved: {file} to {destination_directory}")

    print("All selected files have been moved.")

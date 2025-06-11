#TODO: resets downloaded bed files (deletes them)

import os

def main():
    # delete files in ../data/bed/

    bed_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'bed')
    if os.path.exists(bed_dir):
        for filename in os.listdir(bed_dir):
            file_path = os.path.join(bed_dir, filename)
            if os.path.isfile(file_path):
                os.remove(file_path)

    #delete files in ../data/metadata/

    metadata_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'metadata')
    if os.path.exists(metadata_dir):
        for filename in os.listdir(metadata_dir):
            file_path = os.path.join(metadata_dir, filename)
            if os.path.isfile(file_path):
                os.remove(file_path)

    return


if __name__ == "__main__":
    main()
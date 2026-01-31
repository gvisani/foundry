
import os
import gzip
import shutil

def uncompress(filepath: str) -> None:
    '''
    Uncompresses .gz file

    Deletes compressed file afterwards
    '''
    with gzip.open(filepath, 'rb') as f_in:
        with open(filepath[:-3], 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    
    os.remove(filepath) # delete compressed file


def uncompress_in_dir(dirpath: str) -> None:
    '''
    Recursively because it's fun :)
    '''
    for file_or_dir in os.listdir(dirpath):

        file_or_dir_path = os.path.join(dirpath, file_or_dir)

        if os.path.isdir(file_or_dir_path):
            uncompress_in_dir(file_or_dir_path)
        else:
            if file_or_dir_path.endswith('.gz'):
                uncompress(file_or_dir_path)
        


if __name__ == '__main__':

    uncompress_in_dir('./inference_outputs')





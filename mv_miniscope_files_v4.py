import os
import pandas as pd
import re

def make_dir(path):
    if not os.path.exists(path):
        try:
            os.makedirs(path)
        except OSError:
            raise OSError('mkdir failed for path: ' + path)


date_str = '2020_10_08'
rootdir = 'D:\\Prez\\cheeseboard\\2020-10\\habituation'
rootdir = '/home/przemek/neurodata/cheeseboard/2020-10/habituation'
dated_dir = os.path.join(rootdir, date_str, 'trial')
dest_dir = os.path.join(rootdir, date_str, 'trial')
make_dir(dest_dir)

animal_subdirs = os.listdir(dated_dir)
def scopedir2tuple(subdir):
    subdir = subdir.replace('_',' ')
    digits = re.sub("[HMS]", ' ', subdir)
    return [int(x) for x in digits.strip().split(' ')]

for animal_name in animal_subdirs:
    print(animal_name)
    animal_path = os.path.join(dated_dir, animal_name)
    #rec_info = pd.read_csv(os.path.join(trial_dir, 'settings_and_notes.dat'), sep='\\t')

    animal_subdirs = os.listdir(animal_path)
    session_dirs = [d for d in animal_subdirs if d.startswith('Session')]
    session_no = 1
    if len(session_dirs) > 0:
        session_no = int(session_dirs[-1].replace('Session','')) + 1

    print('Animal subdir: ')
    print(animal_subdirs)
    rec_dirs = [d for d in animal_subdirs if d[0].isdigit()]
    sorted(rec_dirs, key=scopedir2tuple)
    print(rec_dirs)
    for rec_dir in rec_dirs:
        trial_dir = os.path.join(animal_path, rec_dir)
        target_dir = os.path.join(animal_path, 'Session' + str(session_no), rec_dir)
        make_dir(target_dir)

        print('Moving files from ' + rec_dir + ' to: ' + target_dir)
        for filename in os.listdir(trial_dir):
            os.rename(os.path.join(trial_dir, filename), os.path.join(target_dir, filename))
        os.rmdir(trial_dir)

#os.rmdir(dated_dir)


import os
import pandas as pd
import re

def make_dir(path):
    if not os.path.exists(path):
        try:
            os.makedirs(path)
        except OSError:
            raise OSError('mkdir failed for path: ' + path)


date_str = '2019-01-01'
rootdir = 'D:\\Prez\\cheeseboard\\2018-10-habituation'
rootdir = '/mnt/DATA/Prez/cheeseboard/2018-10-learning'
dated_dir = os.path.join(rootdir, date_str, 'caimg')
dest_dir = os.path.join(rootdir, date_str, 'mv_caimg')
make_dir(dest_dir)

miniscope_subdirs = os.listdir(dated_dir)
def scopedir2tuple(subdir):
    digits = re.sub("[HMS]", ' ', subdir)
    return [int(x) for x in digits.strip().split(' ')]

sorted(miniscope_subdirs, key=scopedir2tuple)
for subdir in miniscope_subdirs:
    trial_dir = os.path.join(dated_dir, subdir)
    rec_info = pd.read_csv(os.path.join(trial_dir, 'settings_and_notes.dat'), sep='\\t')
    animal = rec_info['animal'][0]
    animal_dir = os.path.join(dest_dir, animal)
    make_dir(animal_dir)

    session_dirs = os.listdir(animal_dir)
    session_no = 1
    if len(session_dirs) > 0:
        session_no = int(session_dirs[-1].replace('Session','')) + 1

    target_dir = os.path.join(animal_dir, 'Session' + str(session_no), subdir)
    make_dir(target_dir)

    print('Moving files from ' + subdir + ' to: ' + target_dir)
    for filename in os.listdir(trial_dir):
        os.rename(os.path.join(trial_dir, filename), os.path.join(target_dir, filename))
    os.rmdir(trial_dir)

os.rmdir(dated_dir)

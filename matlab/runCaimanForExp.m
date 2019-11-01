expMonth = '2019-08';
expDir = ['~/neurodata/cheeseboard-down/down_2/' expMonth];
animal = 'E-TR';

expTitles = {'habituation', 'learning'};
for expTitle = expTitles
    rootDir = [expDir filesep expTitle{1}];
    subdirs = dir(rootDir);
    for i = 1:numel(subdirs)
        dateStr = subdirs(i).name;
        if startsWith(dateStr, expMonth)
            loadCaimanTrial
        end
    end

end
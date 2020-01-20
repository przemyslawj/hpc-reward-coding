expMonth = '2019-07';
expDir = ['/home/prez/neurodata/cheeseboard-down/down_2/' expMonth];
%expDir = ['/mnt/DATA/Prez/cheeseboard-down/down_2/' expMonth];
animal = 'D-BR';

expTitles = {'habituation', 'learning'};
for expTitle = expTitles
    rootDir = [expDir filesep expTitle{1}];
    subdirs = dir(rootDir);
    for i = 1:numel(subdirs)
        dateStr = subdirs(i).name;
        if startsWith(dateStr, expMonth(1:2))
            loadCaimanTrial
        end
    end

end

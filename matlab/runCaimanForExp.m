%% Merges tracking and caimg traces produced by CaImAn.
% The script requires that the variables: animal and expMonth are present in
% the environment

animal = 'P-BR';
expMonth = '2020-10';
expDir = ['/mnt/DATA/Prez/cheeseboard-down/down_2/' expMonth];
%expDir = ['/mnt/DATA/Prez/cheeseboard/' expMonth];
overwrite = 0;
v3 = 0;

expTitles = {'habituation', 'learning'};
for expTitle = expTitles
    rootDir = [expDir filesep expTitle{1}];
    subdirs = dir(rootDir);
    for i = 1:numel(subdirs)
        dateStr = subdirs(i).name;

        if startsWith(dateStr, expMonth(1:2))
            result_fpath = fullfile(rootDir, dateStr, 'caiman', animal, ...
                                    'filtered', 'traces_and_positions.csv');
			if overwrite || ~isfile(result_fpath)
				loadCaimanTrial
			end
        end
    end

end

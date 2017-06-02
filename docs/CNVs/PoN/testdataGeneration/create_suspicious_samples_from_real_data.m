% THIS FILE IS FOR RECORD KEEPING ONLY
%  THERE IS NO SUPPORT FOR THIS FILE AND IT IS MEANT FOR GATK DEV ONLY
%  For usage of this data, see CoveragePoNQCUtilsUnitTest.java
% LTL January 6, 2016
%
% Please excuse code duplication and other sloppiness.

%% Anonymize the real input data
% This file cannot be made publicly available and will be anonymized.  This
%  input file is not in the GATK CNV repo.
tbl = readtable('no_events_tn.txt', 'FileType', 'text', 'Delimiter', '\t');
num_targets = length(tbl.sample_1);
target_names = cell(num_targets,1);
for i = 1:length(target_names)
    target_names{i} = ['target_' num2str(i)];
end

contig_names = cell(num_targets,1);
for i = 1:length(target_names)
    contig_names{i} = '2';
end

start_pos = round(linspace(1000000, 2000000, num_targets))';
end_pos = start_pos + 50;

% Create a few fake samples that have no events and just a little bit of
%   noise
fake_no_event_samples = log2(2 .^ tbl{:,5:end} + normrnd(0, 0.1, size(tbl{:,5:end})));
added_sample_names = cell(size(fake_no_event_samples,2),1);
for i = 1:length(added_sample_names)
    added_sample_names{i} = ['sample_' num2str(i + width(tbl) - 4)];
end

% Update the table with the anonymized data and new fake data
tbl_fake_data = array2table(fake_no_event_samples, 'VariableNames', added_sample_names);
tbl = [tbl tbl_fake_data];
tbl.NAME = target_names;
tbl.CONTIG = contig_names;
tbl.START = start_pos;
tbl.END = end_pos;

% The original file is anonymized and can be written to disk.  This file is
% in the GATK CNV repo
writetable(tbl,'no_events_tn_an.txt','FileType','text', 'Delimiter','\t');

%% Create artificial amp events in samples 1-3
tbl = readtable('no_events_tn_an.txt', 'FileType', 'text', 'Delimiter', '\t');
tbl.sample_1 = log2((2 .^ tbl.sample_1) + .3);

indicesHalf = 1:floor(length(tbl.sample_2)/2);
tbl.sample_2(indicesHalf) = log2((2 .^ tbl.sample_2(indicesHalf)) + .3);

indicesThird = 1:floor(length(tbl.sample_3)/3);
tbl.sample_3(indicesThird) = log2((2 .^ tbl.sample_3(indicesThird)) + .3);

writetable(tbl,'events_tn.txt','FileType','text', 'Delimiter','\t');

%% Create artificial del events in samples 1-3
tbl2 = readtable('no_events_tn_an.txt', 'FileType', 'text', 'Delimiter', '\t');

tbl2.sample_1 = log2((2 .^ tbl2.sample_1) - .3);

indicesHalf = 1:floor(length(tbl2.sample_2)/2);
tbl2.sample_2(indicesHalf) = log2((2 .^ tbl2.sample_2(indicesHalf)) - .3);

indicesThird = 1:floor(length(tbl2.sample_3)/3);
tbl2.sample_3(indicesThird) = log2((2 .^ tbl2.sample_3(indicesThird)) - .3);

writetable(tbl2,'del_events_tn.txt','FileType','text', 'Delimiter','\t');


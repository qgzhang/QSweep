clear;

dataset{1} = 'alcatraz_courtyard'; %

fp = fopen('spreadsheet.txt','w');
s20 = '--------------------'; 
s8 = '--------'; 
s16 = '----------------'; 
seg0 = sprintf(' %20s-%20s-%20s \n', s20,s20,s20);
seg1 = sprintf('|%20s|%20s|%20s|\n', s20,s20,s20);

blk_size = 100;

files = dir( fullfile('test_results',strcat(dataset{1},'_rec*')) );

time_pl_all = zeros(length(files),1);
mres_pl_all = zeros(length(files),blk_size);
nitr_pl_all = zeros(length(files),1);

time_sw_2v_all = zeros(length(files),1);
mres_sw_2v_all = zeros(length(files),blk_size);
nitr_sw_2v_all = zeros(length(files),1);

t = 0;
XArray_pl_all        = [];
XArray_sw_2v_all     = [];

time_pl = 0;
time_sw_2v = 0;

mres_pl = 0;
mres_sw_2v = 0;

nitr_pl = 0;
nitr_sw_2v = 0;

X_pl = zeros(3,blk_size);
X_sw_2v = zeros(3,blk_size);

for f=1:length(files)

    load(fullfile('test_results',files(f).name));

    time_pl_all(f) = time_pl;
    time_sw_2v_all(f) = time_sw_2v;
   
    t = t + length(mres_pl);
    mres_pl_all(f, 1:length(mres_pl))               = mres_pl;
    mres_sw_2v_all(f, 1:length(mres_sw_2v))         = mres_sw_2v;
   
    nitr_pl_all(f) = nitr_pl;
    nitr_sw_2v_all(f) = nitr_sw_2v;
   
    XArray_pl_all        = [XArray_pl_all        , X_pl];
    XArray_sw_2v_all     = [XArray_sw_2v_all     , X_sw_2v];
    
end

save( strcat('X_and_res_',dataset{1}), ...
        'XArray_pl_all'       , 'mres_pl_all'        , ...
        'XArray_sw_2v_all'    , 'mres_sw_2v_all' );

fprintf(fp, '\n\n');
fprintf(fp,seg0);
fprintf(fp, '%40s\n\n', dataset{1});
fprintf(fp,seg1);
fprintf(fp,'|%20s|%20s|%20s| \n', 'Algorithm', 'Total Runtime', 'Avg Converged Error');
fprintf(fp,seg1);


fprintf(fp,'|%20s|%20.3f|%20.3f|\n', ...
    'Polyhedron Collapse', sum(time_pl_all), sum(sum(mres_pl_all))/t);
fprintf(fp,seg1);

fprintf(fp,'|%20s|%20.3f|%20.3f|\n', ...
    'Q-sweep', sum(time_sw_2v_all), sum(sum(mres_sw_2v_all))/t);
fprintf(fp,seg1);
fclose(fp);

fprintf('\n\n\n');
fprintf('\n\n');
fprintf(seg0);
fprintf('%40s\n\n', dataset{1});
fprintf(seg1);
fprintf('|%20s|%20s|%20s| \n', 'Algorithm', 'Total Runtime', 'Avg Converged Error');
fprintf(seg1);


fprintf('|%20s|%20.3f|%20.3f|\n', ...
    'Polyhedron Collapse', sum(time_pl_all), sum(sum(mres_pl_all))/t);
fprintf(seg1);

fprintf('|%20s|%20.3f|%20.3f|\n', ...
    'Q-sweep', sum(time_sw_2v_all), sum(sum(mres_sw_2v_all))/t);
fprintf(seg1);
fprintf('\n\n\n');


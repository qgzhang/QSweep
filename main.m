
cd('functions'); 

fprintf('\n');
fprintf('\n');

fprintf('******        Demo for ICCV 2017 submission:        *******\n');
fprintf('* Quasiconvex Plane Sweep for Triangulation with Outliers *\n');
fprintf('******    Takes 4 minutes on i7-3.6GHz, 16GB PC:    *******\n');

% fprintf('\n');
% % compile C code if necessary
% try
%     % Load point matches
%     load( 'test_mex.mat' );
%     fprintf('Checking mex ...\n');
%     [~,~] = c_stepsize(aa,bb,cc,dd,X,dir,K);
%     fprintf('Mex check passed.\n');
% catch exception
%     fprintf('Compiling source code...\n');
%     mex c_stepsize.c
% end




clear;
dataset{1} = 'alcatraz_courtyard'; 
data_path = fullfile('dataset', dataset{1}, 'data');
load(data_path);

M = size(PArray_r,1)/3;
N = size(ImgArray_r,2);

P_all = PArray_r;
Img_all = ImgArray_r;
Norm = 3;

bound = N;
blk_size = 100;

fprintf('\n');
fprintf('dataset : %s \n','Alcatraz Courtyard');
fprintf('\n');

fprintf('running...\n');

for blk = 1:1000
    
    pts = blk_size*(blk-1)+1 : blk_size*blk;
    pts(pts>bound) = [];

    if isempty(pts)
        break;
    end
    save_path = fullfile('test_results', strcat(dataset{1},...
        '_rec_sw_', num2str(pts(1)) ,'_' , num2str(pts(end))));

    
    %%%%%  Polyhedron collapse %%%%%%%%%%%%%%%%%%%%
    tic
    [ X_pl, mres_pl, nitr_pl] = f_poly(P_all, Img_all(:,pts));
    time_pl = toc;

%plot3(X_pl(1,:),X_pl(2,:),X_pl(3,:),'.');
    %%%%%  Q-sweep %%%%%%%%%%%%%%%%%%%%
    tic 
    % initialise X by the 2-views method
    [ X_2v ] = f_2views_fea( P_all, Img_all(:,pts) );
    [ X_sw_2v, mres_sw_2v, nitr_sw_2v] = ...
        f_Qsweep(P_all, Img_all(:,pts), X_2v );    
    time_sw_2v = toc;

    save(save_path, ...
        'X_pl',        'mres_pl',        'nitr_pl',        'time_pl', ...
        'X_sw_2v',     'mres_sw_2v',     'nitr_sw_2v',     'time_sw_2v');

    
    if mod(blk,10)==0
        fprintf('    %.2f%% completed\n',100*blk/295);
    end
    
end

%%%%%  summarise converged error and run time %%%%%%%%%%%%%%%%%%%%
res_and_time;

%%%%%  compare 3D structures %%%%%
fig_3D;

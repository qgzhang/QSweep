clear;

dataset{1} = 'alcatraz_courtyard'; %

eps = 1e-6;
varname{1} = 'pl';
varname{2} = 'sw_2v';

dataname{1} = 'Courtyard'; %

colors{1} = 'yellow';
colors{2} = 'red';

axis_config{1} = [-0.15, 0.3, -0.17, 0.17, 0, 0.6];

position{1} = [0,0, 900, 900, 0, 90];

bin{1} = 8;

load( strcat('X_and_res_',dataset{1}) );


%% 3D structure by ploy  
markersize{1} = 8;

X = XArray_pl_all;
id1 = X(1,:) >= axis_config{1}(1) & X(1,:) <= axis_config{1}(2);
id2 = X(2,:) >= axis_config{1}(3) & X(2,:) <= axis_config{1}(4);
id3 = X(3,:) >= axis_config{1}(5) & X(3,:) <= axis_config{1}(6);
id_pl = id1 & id2 & id3;

X = XArray_sw_2v_all;
id1 = X(1,:) >= axis_config{1}(1) & X(1,:) <= axis_config{1}(2);
id2 = X(2,:) >= axis_config{1}(3) & X(2,:) <= axis_config{1}(4);
id3 = X(3,:) >= axis_config{1}(5) & X(3,:) <= axis_config{1}(6);
id_sw = id1 & id2 & id3;

f1 = figure; 
X = XArray_pl_all;
X = X(:,id_pl(1:size(X,2)));
plot3(X(1,:),X(2,:),X(3,:),'b.','markersize',markersize{1}); hold on;
%     pcx = pointCloud(X');
%     pcsetcolor(pcx, 'b');
%     pcshow(pcx,'markersize',40);
%     xlabel('x');
%     ylabel('y');
%     zlabel('z');
axis(axis_config{1});
view(position{1}(5:6));
set(f1, 'Position', position{1}(1:4));
%     set(gca,'visible','off');
%     set(gca, 'XTICK', []);  
%     set(gca, 'YTICK', []);
%     set(gca, 'ZTICK', []);
grid on;
set(gca,'XTickLabel','','YTickLabel','','ZTickLabel','');
title(strcat(dataname{1},' - Polyhedron Collapse'),'fontsize',18);
% savefig(f1, fullfile('figures', strcat('reconst_Poly_',dataname{1}) ));

f2 = figure; 
X = XArray_sw_2v_all;
X = X(:,id_sw(1:size(X,2)));
plot3(X(1,:),X(2,:),X(3,:),'b.','markersize',markersize{1}); hold on;
%     pcx = pointCloud(X');
%     pcsetcolor(pcx, 'b');
%     pcshow(pcx,'markersize',40);
%     xlabel('x');
%     ylabel('y');
%     zlabel('z');
axis(axis_config{1});
view(position{1}(5:6));
set(f2, 'Position', position{1}(1:4));
%     set(gca,'visible','off');
title(strcat(dataname{1},' - Q-sweep'), 'fontsize',18);
%     set(gca, 'XTICK', []);
%     set(gca, 'YTICK', []);
%     set(gca, 'ZTICK', []);
grid on;
set(gca,'XTickLabel','','YTickLabel','','ZTickLabel','');
% savefig(f2, fullfile('figures', strcat('reconst_Sweep_',dataname{1}) ));
clear;





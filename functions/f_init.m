function [ XArray ] = f_init( PArray, ImgArray)

% initialise X to be infront of all visible cameras
% rmpath('C:\Program Files\Mosek\8\toolbox\r2014a');
% rmpath( fullfile('..', '..', '..', 'libraries', 'mosek', '8', 'toolbox', 'r2014a') );
N = size(ImgArray,2);
M = size(PArray,1)/3;

eps = 1e-7;
XArray = zeros(3,N);
[A1,A2,B1,B2,C,D] = f_gen_coef_2(PArray, ImgArray);

options = optimoptions('linprog','Display','off');

for n = 1:N

    a1 = A1(M*(n-1)+1:M*n,:);
%     a2 = A2(M*(n-1)+1:M*n,:);
%     b1 = B1(M*(n-1)+1:M*n,:);
%     b2 = B2(M*(n-1)+1:M*n,:);
    c  = C;
    d  = D;
    
    fid = find(isfinite(sum(a1,2)));
%     a1 = a1(fid,:);
%     a2 = a2(fid,:);
    c  = c (fid,:);
%     b1 = b1(fid,:);
%     b2 = b2(fid,:);
    d  = d (fid,:);
    
%     aa = [a1;-a1;a2;-a2];
%     bb = [b1;-b1;b2;-b2];
%     cc = [c;c;c;c];
%     dd = [d;d;d;d];
    
    AA = -c;
    BB = d-eps;
    X_n = linprog([0;0;0], AA,BB, [],[],[],[],[], options);
    
    XArray(:,n) = X_n;    

end % end of for


end


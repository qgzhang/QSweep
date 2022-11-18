function [ XArray, medianres, nitr] = f_Qsweep(PArray, ImgArray, Xest )

% compare the ransac and the sweeping initialised with ransac

N = size(ImgArray,2);
XArray = zeros(3,N);

M = size(PArray,1)/3;

medianres = zeros(1,N);

[A1,A2,B1,B2,C,D] = f_gen_coef_2(PArray, ImgArray);
nitr = 0;

for n = 1:N
%     fprintf('sweep, n=%d\n',n);
    a1 = A1(M*(n-1)+1:M*n,:);
    a2 = A2(M*(n-1)+1:M*n,:);
    b1 = B1(M*(n-1)+1:M*n,:);
    b2 = B2(M*(n-1)+1:M*n,:);
    c  = C;
    d  = D;
    
    fid = find(isfinite(sum(a1,2)));
    a1 = a1(fid,:);
    a2 = a2(fid,:);
    c  = c (fid,:);
    b1 = b1(fid,:);
    b2 = b2(fid,:);
    d  = d (fid,:);
    
    aa = [a1;-a1;a2;-a2];
    bb = [b1;-b1;b2;-b2];
    cc = [c;c;c;c];
    dd = [d;d;d;d];
    
    X0 = Xest(:,n);
    
    [X, mres_s, nitr_s] = f_sweep_single(aa,bb,cc,dd,X0);
 
    nitr = nitr + nitr_s;
    XArray(:,n) = X;
    medianres(n) = mres_s;

end


end

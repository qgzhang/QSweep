function [ XArray, medianres, nitr] = f_poly(PArray, ImgArray)

N = size(ImgArray,2);
XArray = zeros(3,N);
medianres = zeros(1,N);
M = size(PArray,1)/3;

[A1,A2,B1,B2,C,D] = f_gen_coef_2(PArray, ImgArray);
nitr = 0;

for n = 1:N
    
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

    id1 = (fid-1)*3+1;
    id2 = (fid-1)*3+2;
    id3 = (fid-1)*3+3;
    pid = [id1';id2';id3'];
    pid = pid(:);
    id1 = (fid-1)*2+1;
    id2 = (fid-1)*2+2;
    iid = [id1';id2'];
    iid = iid(:);

    
    X = f_2views_fea_single( PArray(pid,:), ImgArray(iid,n), c, d );
    [X, mres_s, nitr_s] = f_poly_single(aa, bb, cc, dd, X);
    
    
%     tic
%     X2 = f_tri_svd( PArray(pid,:), ImgArray(iid,n) );
%     toc
%     
%     [X,X2]
%     
%     res1 = (f_res(X,PArray(pid,:), ImgArray(iid,n),3));
%     res2 = (f_res(X2,PArray(pid,:), ImgArray(iid,n),3));
%     
%     [res1, res2]
    
    nitr = nitr + nitr_s;
    XArray(:,n) = X;
    medianres(n) = mres_s;

end % for n=1:N

end





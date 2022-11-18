function [ roots ] = f_root( aa, bb, cc, dd, i, j, X, dir, lb, ub )

%     ma = aa(i,:);
%     mb = bb(i,:);
%     mc = cc(i,:);
%     md = dd(i,:);
% 
%     na = aa(j,:);
%     nb = bb(j);
%     nc = cc(j,:);
%     nd = dd(j);
% 
%     a = nc*(ma*dir*dir) - na*(dir*(mc*dir));
%     b = nc*(ma*X*dir) + nc*(ma*dir*X) + nd*(ma*dir) + nc*(mb*dir) ...
%             -na*(X*mc*dir) - na*(dir*mc*X) - na*(md*dir) - nb*(mc*dir);
%     c = nc*(ma*X*X) + nd*(ma*X) + nc*(mb*X) + nd*md ...
%         -na*(X*mc*X) - na*(md*X) - nb*(mc*X) - nb*md;

    eps = 1e-15;
    
    a1 = aa(i,:);
    b1 = bb(i,:);
    c1 = cc(i,:);
    d1 = dd(i,:);

    a2 = aa(j,:);
    b2 = bb(j);
    c2 = cc(j,:);
    d2 = dd(j);
    
    

% fea = (cc([i,j],:)*X+dd([i,j]))   
% res = (aa([i,j],:)*X+bb([i,j])) ./ (cc([i,j],:)*X+dd([i,j]))
% X1 = X + 1e-10*dir;
% res1 = (aa([i,j],:)*X1+bb([i,j])) ./ (cc([i,j],:)*X1+dd([i,j]));
% incremental = res1-res


    a = (a1*dir)*(c2*dir) - (a2*dir)*(c1*dir);
    b = (a1*X)*(c2*dir) + (a1*dir)*(c2*X) + (a1*dir)*d2 + b1*(c2*dir) ...
       -(a2*X)*(c1*dir) - (a2*dir)*(c1*X) - (a2*dir)*d1 - b2*(c1*dir);
    c = (a1*X)*(c2*X) + (a1*X)*d2 + b1*(c2*X) + b1*d2 ...
       -(a2*X)*(c1*X) - (a2*X)*d1 - b2*(c1*X) - b2*d1;

    roots = [(-b+sqrt(b.*b-4*a.*c))./(2*a) , (-b-sqrt(b.*b-4*a.*c))./(2*a)];

    
    
    %fea_1 = -(c1*X + d1) / (c1*dir);
    %fea_2 = -(c2*X + d2) ./ (c2*dir);
    %Fea_2 = repmat(fea_2,2,1);
    
%     id = find(j==i);
%     roots(id,:) = [-1,-1];
    roots = roots(:);
    
    
    for i=1:length(roots) % eliminate imaginary roots 
        if ~isreal(roots(i))
            roots(i) = -1;
        end
    end
    roots(~isfinite(roots)) = -1;% roots with i itself is INF
    
    
    % delete infeasible
    if exist('lb','var') && exist('ub','var')
        roots(roots<lb) = -1;
        roots(roots>ub) = -1;
    else
        roots( c1*dir*roots + c1*X + d1 <= eps ) = -1; 
        roots( repmat(c2,2,1) * dir .* roots  +  repmat(c2,2,1) * X + repmat(d2,2,1) <= eps ) = -1;
    end
    
    % delete negative intersections
    
%     values = (a1*X+a1*dir*roots + b1) ./ (c1*X+c1*dir*roots + d1);
%     neg = (values<=eps);
%     roots(neg) = -1;
    
    
end


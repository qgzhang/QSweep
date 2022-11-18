function [ X, mres, nitr] = f_poly_single(aa, bb, cc, dd, X)


    mres0 = Inf;
    mres = 1e10;
    flag = true;

    nitr = 0;

    while abs(mres0-mres)>1e-16
        nitr = nitr+1;
        % find the active constraints

        res = (aa * X + bb) ./ (cc * X + dd);
        %tt = [res(1:M),res(M+1:2*M),res(2*M+1:3*M),res(3*M+1:4*M)];
        mres0 = mres;
        mres = max(res);
        if mres>mres0
            flag = false;
            break;
        end
        Active = find(mres-res<=1e-4*mres)';

        % improving direction
        %  - the normal of a plane ax+by+cz+d = 0 is [a,b,c]
        if length(Active)==1
            n1 = -( (aa(Active(1),:) - mres*cc(Active(1),:)) );
            n1 = n1./norm(n1);
            dir = n1;
        elseif length(Active)==2
            n1 = -( (aa(Active(1),:) - mres*cc(Active(1),:)) );
            n2 = -( (aa(Active(2),:) - mres*cc(Active(2),:)) );
            n1 = n1./norm(n1);
            n2 = n2./norm(n2);
            dir = n1+n2;

        elseif length(Active)==3
            n1 = -( (aa(Active(1),:) - mres*cc(Active(1),:)) );
            n2 = -( (aa(Active(2),:) - mres*cc(Active(2),:)) );
            n3 = -( (aa(Active(3),:) - mres*cc(Active(3),:)) );
            n1 = n1./norm(n1);
            n2 = n2./norm(n2);
            n3 = n3./norm(n3);
            dp = cross(n1,n2) + cross(n2,n3) + cross(n3,n1);
            % check coplanar
            dp = f_coplanar(dp,n1,n2,n3);
            if nnz(dp)==0
                break;
            end
            s = sign(dot(n1,dp)); % unsure about dp
            dir = s*dp;


        elseif length(Active)>3
            triplet = nchoosek(Active,3);
            for i = 1:size(triplet,1)
                tri = triplet(i,:);
                n1 = -(aa(tri(1),:) - mres*cc(tri(1),:));
                n2 = -(aa(tri(2),:) - mres*cc(tri(2),:));
                n3 = -(aa(tri(3),:) - mres*cc(tri(3),:));
                n1 = n1./norm(n1);
                n2 = n2./norm(n2);
                n3 = n3./norm(n3);
                dp = cross(n1,n2) + cross(n2,n3) + cross(n3,n1);
                dp = f_coplanar(dp,n1,n2,n3);
                if nnz(dp)==0
                    continue;
                end
                s = sign(dot(n1,dp)); % unsure about dp
                dir = s*dp;
                valid = true;
                for j = setdiff(Active,tri)
                    nj = -(aa(j,:) - mres*cc(j,:));
                    if sign(dot(dir,nj)) == -1 % oppsite direction: constraint j will increase along dir
                        valid = false;
                        break; %for j = setdiff(Active,tri)  ---not a good dir
                    end
                end
                if valid==true
                    % found a valid direction
                    break; %for i = size(triplet,1)
                end
            end
            if valid == false
                % non valid direction, so terminate triangulating the point n
                break; % while
            end
        else
            break;
        end % end of if
        dir = dir'/norm(dir);  % as a unit col vector

        %% step size


        %method = 1; % least positive root
        method = 2; % best reduction root
        %method = 3; % optimal stepsize
        alpha1 = stepsize(Active,method);
        if isempty(alpha1) || alpha1==0 % no positive roots implies either convergence or the direction is incorrect
            %dir = redir();
            break;
        end


        %% update X
        alpha = alpha1;
        X0 = X;
        %max(f_res(X, PArray, ImgArray,3))
        X = X + alpha * dir;
        %max(f_res(X, PArray, ImgArray,3))

    end % while

    if flag
        X = X;
    else
        X = X0;
    end

    nCam = length(dd)/4;
    K = ceil(nCam/2); 
    res_cons = (aa*X+bb)./(cc*X+dd);
    res_cam = reshape(res_cons,nCam,4);
    res_cam = max(res_cam,[],2);
    [sres,~] = sort(res_cam,'descend');
    %mres = sres(K); % median error
    mres = max(res_cons); % max error (converged error)








    function [dir,flag] = redir(dir)
        flag = 0; % unable to update dir
        if length(Active)==1
            return;
        elseif length(Active)==2
            if sign(dot(n1,n2)) == 1 % both n1 and n2 are valid 
                flag = 1;
                dir = n1;
            else
                return;
            end
        elseif length(Active)==3
            if sign(dot(n1+n2,n3)) == 1 % valid
                flag = 1;
                dir = n1+n2;
            elseif sign(dot(n2+n3,n1)) == 1 % valid
                flag = 1;
                dir = n2+n3;
            elseif sign(dot(n3+n1,n2)) == 1 % valid
                flag = 1;
                dir = n3+n1;
            else
                return;
            end
        elseif length(Active)>3
            triplet = nchoosek(Active,3);
            triplet(1:end,:) = triplet(end:-1:1,:);
            for k = 1:size(triplet,1)
                tri = triplet(k,:);
                n1 = -(aa(tri(1),:) - mres*cc(tri(1),:));
                n2 = -(aa(tri(2),:) - mres*cc(tri(2),:));
                n3 = -(aa(tri(3),:) - mres*cc(tri(3),:));
                n1 = n1./norm(n1);
                n2 = n2./norm(n2);
                n3 = n3./norm(n3);
                dp = cross(n1,n2) + cross(n2,n3) + cross(n3,n1);
                dp = f_coplanar(dp,n1,n2,n3);
                if nnz(dp)==0
                    continue;
                end
                s = sign(dot(n1,dp)); % unsure about dp
                dp = s*dp;
                valid = true;
                for l = setdiff(Active,tri)
                    nl = -(aa(l,:) - mres*cc(l,:));
                    if sign(dot(dir,nl)) == -1 % oppsite direction: constraint j will increase along dir
                        valid = false;
                        break; %for j = setdiff(Active,tri)  ---not a good dir
                    end
                end
                if valid==true
                    % found a valid direction
                    break; %for i = size(triplet,1)
                end
            end
            if nnz(dp-dir<1e-6)==3 
                return
            else
                flag=1;
                dir = dp;
            end
        end
    end % end of function redir()

    function step = stepsize(Active, method)
        if ~ismember(method,[1,2,3])
            printf('step size method is not recognised\n');
            step = [];
        end
        if method == 3
            % method 3: get optimal stepsize by bisection
            ua=10000000;
            la=0;
            step = 0;
            while (ua-la)>1e-7
                mida = (ua+la)/2;
                X3_m = X + (mida) * dir;
                X3_l = X + (mida-1e-6) * dir;
                X3_u = X + (mida+1e-6) * dir;
                if sum(cc * X3_m + dd < 1e-8)>0
                    ua = mida;
                    continue;
                end
                mres_m = max((aa * X3_m + bb) ./ (cc * X3_m + dd));
                mres_l = max((aa * X3_l + bb) ./ (cc * X3_l + dd));
                mres_u = max((aa * X3_u + bb) ./ (cc * X3_u + dd));
                if (mres_u > mres_m) && (mres_m > mres_l) 
                    ua = mida;
                elseif (mres_l > mres_m) && (mres_m > mres_u) 
                    la = mida;
                else
                    break;
                end
            end
            step = mida;
        else
            a0 = aa(Active,:);
            b0 = bb(Active,:);
            c0 = cc(Active,:);
            d0 = dd(Active,:);
            act_dev= ( (c0*X+d0).*(a0*dir) - (a0*X+b0).*(c0*dir) ) ./ ((c0*X+d0).^2);
            % master constraint: the least negative one, so is the maximum one
            [~,ind] = max(act_dev);
            master = Active(ind);
            Roots = f_root(aa,bb,cc,dd,master,1:size(aa,1),X,dir);
            if method == 1
                step = min(Roots(Roots>1e-10)); 
            elseif method == 2
                minmax = mres;
                brt = 0;
                rts = Roots(Roots>0);
                for r = 1:length(rts)
                    Xp = X + rts(r)*dir;
                    resp = (aa*Xp+bb)./(cc*Xp+dd);
                    if minmax > max(resp);
                        minmax = max(resp);
                        brt = rts(r);
                    end
                end
                step = brt;
            end
        end
    end


end




function dp = f_coplanar(dp,n1,n2,n3)         
% test combination of each pair normals
    if( abs(dot(n1,dp))<1e-6 && abs(dot(n2,dp))<1e-6 && abs(dot(n3,dp))<1e-6 )
        if sign(dot(n1+n2,n3)) == 1 % valid
            dp = n1+n2;
        elseif sign(dot(n2+n3,n1)) == 1 % valid
            dp = n2+n3;
        elseif sign(dot(n3+n1,n2)) == 1 % valid
            dp = n3+n1;
        else
            dp = zeros(size(dp));
        end
    end
end



function [ X, mres_s, nitr] = f_sweep_single(aa, bb, cc, dd, Xest)

    X = Xest;
  
    nCam = length(dd)/4;
    K = ceil(nCam/2); % the median camera

    mres0 = Inf;
    mres = 1e10;
    flag = true;
    nitr = 0;

    while abs(mres0-mres)>1e-6
        nitr = nitr+1;
%         fprintf('nitr = %d\n',nitr);
        % find the active constraint: 
        % the constraint with error equaling to the error of median camera
        res_cons = (aa * X + bb) ./ (cc * X + dd);
        res_cam = reshape(res_cons,nCam,4);
        res_cam = max(res_cam,[],2);
        mres0 = mres;
        [sres,~] = sort(res_cam,'descend');
        mres = sres(K); % the median camera's error
        if mres>mres0
            flag = false;
            break;
        end
        %Active = find(abs(mres-res)<=eps)';
        Active = find(abs(mres-res_cons)<=1e-4*mres)';
        
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
                    if sign(dot(dir,nj)) == -1
                        valid = false;
                        break; %for j = setdiff(Active,tri) 
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
        dir = dir'/norm(dir)/100;  % divided by 100 to increase numerical stability
        
        %% step size
        
        [alpha,~] = c_stepsize(aa, bb, cc, dd, X, dir, K);
        if alpha <= 1e-8
            break;
        end
        
        X0 = X;
        X = X + alpha * dir;
    end % while
    
    if flag == false
        X = X0;
        mres_s = mres0;
    else
        mres_s = mres;
    end
    
end % end of function f_poly_LMS_C


% check coplanarity
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

function showGeom(no,el,nopr)
    
    %to do:
    %return handles (optionally) to user?
    %plot applied moments
    %take force and moment vectors into account for deciding tag placement
    %process special nodes
    %(z-index visibility thing for annotations)

    %do not allow running this script directly; instead, call with input
    %from a script
    if nargin == 0; error('Do not run directly; call from a script instead.'); end

    %calculate bounding box based on all visual elements (labels included)
    %get normal vector (and catch error) in case of 2d mech
    %visualize boundary conditions
    %visualize loads and inputs
    
    figure;
    
    %legend axes
    subplot(1,2,1);
    hold on
    leg1 = plot(0,0); text(1,0,'node');
    leg2 = plot(0,-0.2); text(1,-0.2,'node with mass and/or inertia');
    leg3 = plot(0,-0.4); text(1,-0.4,'fixed node (all directions)');
    leg4 = plot(0,-0.6); text(1,-0.6,'fixed node (one or two dir.)');
    leg5 = text(0,0.2,'#'); text(1,0.2,'node number');
    leg6 = text(0,0.4,'#'); text(1,0.4,'element number');

    set(gca,'DataAspectRatio',[5 1 1])
    set(gca,'Position',[0.1 0.1 0.3 0.8])
    set(gca,'Visible','off')
    set(gca,'XLim',[0 2])
    set(gca,'YLim',[-1 1])
    
    %mechanism axes
    ah = subplot(1,2,2);
    hold on
    
    %handy vars
    nno = size(no,1);
    nel = size(el,1);
    
    %check input
        % no input and load on same node-direction
        
        % elements are connected to defined nodes
        if any(el(:)<=0), error('Element seems connected to node number <=0.'); end
        if max(el(:))>nno, error('Element seems connected to node that does not exist.'); end
        if any((el(:,1)-el(:,2))==0), error('Element seems connected to the same node.'); end
        
    %plot elements
    eh = [];
    for i=1:nel
        xp = no(el(i,1),:);
        xq = no(el(i,2),:);

%         %check element length: if not too short, leave some space for the
%         %element number tag
%         minlen = ;
%         d = 0.1;
%         len = norm(xq - xp);
%         if len > minlen
%             x1 = xp+1/len*(xq-xp)*(len/2-d/2);
%             x2 = xp+1/len*(xq-xp)*(len/2+d/2);
%             eh(end+1) = plot3(...
%                 [xp(1) x1(1)],[xp(2) x1(2)],[xp(3) x1(3)]);
%             eh(end+1) = plot3(...
%                 [x2(1) xq(1)],[x2(2) xq(2)],[x2(3) xq(3)]);
%                 
%         else
%             %use index end+1 instead of i, because elements have one or two
%             %line handles, depending on element length
%                 eh(end+1) = plot3(...
%                     [xp(1) xq(1)],[xp(2) xq(2)],[xp(3) xq(3)]);
%         end

          eh(i) = plot3([xp(1) xq(1)],[xp(2) xq(2)],[xp(3) xq(3)]);
    end
    
    %plot nodes
    nh = [];
    for i=1:nno
        nh(i) = plot3(no(i,1),no(i,2),no(i,3));
    end
    
    %get largest and characteristisch distance of system
    maxdistiter = [];
    maxdistiteri = [];
    for i=1:nno
       relvec = no-repmat(no(i,:),size(no,1),1);
       dist = sqrt(sum(relvec.*relvec,2));
       [maxdistiter(i), maxdistiteri(i)] = max(dist);
    end
    [maxdist, maxdisti] = max(maxdistiter);
    nmax2 = maxdistiteri(maxdisti);
    nmax1 = maxdisti;
    clen = maxdist/20; %for distance between node and its label
    arrowlength = maxdist/6; %for length of arrows

    %get average plane through nodes
    [normalvec, mech_dim, normalvec2] = getPlane(no);
        
    %plot node numbers
    nnh = [];
    for i=1:nno
        
        %get unit vectors pointing out of node
        outvectors = [];
        [ni,nj]=find(el==i); %indices of node in el matr
        node_mult = length(ni); %nr of elements attached to node
        
        if node_mult == 0
            npos = no(i,:) + clen*[1 1 1];
        else
            for j=1:node_mult
                if nj(j)==1, k=2; end
                if nj(j)==2, k=1; end
                connecting_node = el(ni(j),k);
                connecting_vec = no(connecting_node,:) - no(i,:);
                outvectors(j,:) = normVec(connecting_vec);
            end
        end
        
        if node_mult == 1
            dir = cross(outvectors(1,:),normalvec);
            npos = no(i,:) + clen*normVec(dir);
        end
        if node_mult == 2 
        
            ang = acosd(dot(outvectors(1,:),outvectors(2,:)));
            if ang < 2
                % parallel
                dir = cross(outvectors(1,:),normalvec);
                npos = no(i,:) + clen*normVec(dir);
            elseif ang > 178
                % anti-parallel
                dir = cross(outvectors(1,:),normalvec);
                npos = no(i,:) + clen*normVec(dir);
            else
                %  span plane, traverse midangle
                bisec = normVec(outvectors(1,:)+outvectors(2,:));
                npos = -clen*bisec + no(i,:);
            end
            
        end
        
        if node_mult == 3
            % bisector of two vecs, then another bisector with remain. 
            bisec1 = normVec(outvectors(1,:)+outvectors(2,:));
            bisec2 = normVec(bisec1+outvectors(3,:));
            npos = -clen*bisec2 + no(i,:);
        end
        
        if node_mult > 3
            %  give up on heuristics, just choose something
            npos = no(i,:) + clen*[1 1 1];
        end
        
        nnh(i) = text(npos(1),npos(2),npos(3),num2str(i));
    end
    
    %plot element numbers
    enh = [];
    for i=1:nel
        enh(i) = text(...
                    [no(el(i,1),1)+no(el(i,2),1)]/2, ...
                    [no(el(i,1),2)+no(el(i,2),2)]/2, ...
                    [no(el(i,1),3)+no(el(i,2),3)]/2, ...
                    num2str(i));
    end
    
    %visual properties
    set([nh leg1],'Marker','.','MarkerSize',16,'Color','k')
    set([leg2],'Marker','.','MarkerSize',16,'Color','r')
    set([leg3],'Marker','.','MarkerSize',16,'Color','b')
    set([leg4],'Marker','.','MarkerSize',16,'Color','g')
    set(eh,'Color',0.75*[1 1 1],'LineWidth',2.5)
    set([nnh leg5],'BackgroundColor','w','Color','b')
    set([enh leg6],'BackgroundColor','w','EdgeColor','k','FontWeight','bold')
    grid on
    
    %axes
    xlabel('x');
    ylabel('y');
    zlabel('z');
    axis(getAxesLim(no,arrowlength))
    set(ah,'DataAspectRatio',[1 1 1])
    if mech_dim == 3
        set(ah,'Projection','perspective')
    end
    
    %set view
    try
        if mech_dim == 1
            [az,el] = view(normalvec);
            [az2,el2] = view(normalvec2);
            azs = [az az2];
            els = [el el2];
            ii = [el el2]>=0;
            switch sum(ii)
                case 0
                    view(az,el)
                case 1
                    view(azs(ii),els(ii))
                case 2
                    [~,mi] = min(els);
                    view(azs(mi),els(mi));
            end
        else
            [az,el] = view(normalvec);
            %do not view mechanism from below horizon, just override
            if el<3
                %TO DO: for 2-D mechanisms, this should be different
                view(az,15)
            end
        end
        
    catch    
        disp('Issue setting viewport.')
    end
        
    %plot force
    for i=1:size(nopr,2)
        if(isfield(nopr(i),'force') && ~isempty(nopr(i).force))
            arrow3(no(i,:)-normVec(nopr(i).force)*arrowlength, ...
                    no(i,:)-normVec(nopr(i).force)*clen/2,[],1.5,3)
        end
    end
    
end

function [normalvec, mech_dim, normalvec2] = getPlane(no)

% OLD CODE:
%     %%%warning off MATLAB:rankDeficientMatrix %surpress rank deficiency warning (happens for 2-D mechanisms)
%     nno = size(no,1);
%     
%     plane_Amat = [no(:,1) no(:,2) ones(nno,1)];
%     plane_Bmat = no(:,3);
%     rank(plane_Amat);
%     if rank(plane_Amat) == 1
%         
%     elseif rank(plane_Amat) == 2
%         error('problem visualizing 2-d mechanisms. Define a normal vector!');
%         return;
%         
%     elseif rank(plane_Amat) == 3
%         planecoef = plane_Amat\plane_Bmat;
%         normalvec = planecoef;
%         normalvec(3) = -1;
%         normalvec = normVec(normalvec);
%     end
    
    %%% NEW CODE:
    %%% from: http://www.ilikebigbits.com/blog/2015/3/2/plane-from-points
    
    %TO DO: check - at least three points required
    %points relative to centroid (average value)
    no_c = no - mean(no,1);
    
    %covariance components (excluding symmetries)
    xx = sum(no_c(:,1).^2);
    yy = sum(no_c(:,2).^2);
    zz = sum(no_c(:,3).^2);
    xy = sum(no_c(:,1).*no_c(:,2));
    xz = sum(no_c(:,1).*no_c(:,3));
    yz = sum(no_c(:,2).*no_c(:,3));
    
    %determinants
    det_x = yy*zz - yz^2;
    det_y = xx*zz - xz^2;
    det_z = xx*yy - xy^2;
    
    det_max = max([det_x det_y det_z]);
        
    %number of zero determinants
    normalvec2 = [];
    switch sum([det_x det_y det_z]==0) 
        case 3
            %nodes lie on a line (1-D mechanism)
            %TO DO: what's a reasonable normalvec here?
            
            %get a vector in vectorspace representing nodes
            v1 = no(2,:) - no(1,:);
            kern = null(v1);
            
            normalvec = kern(:,1);
            normalvec2 = kern(:,2);
            mech_dim = 1;
            return;
        case 2
            %nodes lie on a plane (2-D mechanism)
            mech_dim = 2;
        otherwise
            mech_dim = 3;
    end
    
    %pick path with best conditioning
    switch det_max
        case det_x
            a = (xz*yz - xy*zz)/det_x;
            b = (xy*yz - xz*yy)/det_x;
            ding = [1,a,b];
        case det_y
            a = (yz*xz - xy*zz)/det_y;
            b = (xy*xz - yz*xx)/det_y;
            ding = [a,1,b];
        case det_z
            a = (yz*xy - xz*yy)/det_z;
            b = (xz*xy - yz*xx)/det_z;
            ding = [a,b,1];
    end
    
    normalvec = normVec(ding);
    
end

function res = normVec(vec)
    res = vec/norm(vec);
end

function lim = getAxesLim(no,arrowlength)
    offset = 0.1;
    xmax = max(no(:,1));
    xmin = min(no(:,1));
    dx = xmax - xmin;
    
    ymax = max(no(:,2));
    ymin = min(no(:,2));
    dy = ymax - ymin;
    
    zmax = max(no(:,3));
    zmin = min(no(:,3));
    dz = zmax - zmin;
    
    d = max([dx dy dz]);
    
    corrx = max([arrowlength offset*d]);
    corry = max([arrowlength offset*d]);
    corrz = max([arrowlength offset*d]);
    lim = [xmin xmax ymin ymax zmin zmax] + [-corrx corrx -corry corry -corrz corrz];
end
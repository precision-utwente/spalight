function showGeom(no,el,nprops,eprops)
    
    %%%%%%%%%%
    %visualize user-specified geometry based on the SPACAR Light input data
    %created by Marijn Nijenhuis
    %%%%%%%%%%
    
    %contributions from Wouter van Dijk
    
    %%%
    %to do:
    %%%
    %return handles (optionally) to user?
    %plot applied moments
    %take force and moment vectors into account for deciding tag placement
    %process special nodes
    %(z-index visibility thing for annotations)
    %calculate bounding box based on all visual elements (labels included)
    %get normal vector (and catch error) in case of 2d mech
    %visualize boundary conditions
    %visualize loads and inputs
    
    %do not allow running this script directly; instead, call with input
    %from a script
    if nargin == 0; error('Do not run directly; call from a script instead.'); end
    
    fh = figure;
    
    %--show legend flag--
    showlegend = false;
    
    if showlegend
        %axes for the legend
        lah = subplot(1,2,1);
        hold on

        dx = 0.2;
        dy = 0.2;
        leg6 = text(0,2*dy,'#'); text(dx,2*dy,'element number');
        leg5 = text(0,1*dy,'#'); text(dx,1*dy,'node number');
        leg1 = plot(0,0*dy); text(dx,0*dy,'node');
        leg2 = plot(0,-1*dy); text(dx,-1*dy,'node with mass and/or inertia');
        leg3 = plot(0,-2*dy); text(dx,-2*dy,'fixed node (all directions)');
        leg4 = plot(0,-3*dy); text(dx,-3*dy,'fixed node (one or two dir.)');

        set(lah,'PlotBoxAspectRatio',[1.8 1 1])
        set(lah,'XLim',[-0.2 1.2])
        set(lah,'YLim',[-0.7 0.5])
        set(lah,'xtick',[],'ytick',[])
        box(lah,'on')

        %axes for the mechanism
        ah = subplot(1,2,2);
        hold on
    else
        ah = axes;
        hold on
    end
    
    %rotation object with callback to update label position after rotation
    rotobj = rotate3d;
    rotobj.ActionPostCallback = @rotationcallback;
    rotobj.Enable = 'on';
    
    %handy vars
    nno = size(no,1);   % Number of nodes
    nel = size(el,1);   % Number of elements
    
    %--check input--
    % no input and load on same node-direction
    
    % elements are connected to defined nodes
    if any(el(:)<=0), error('Element seems connected to node number <=0.'); end
    if max(el(:))>nno, error('Element seems connected to node that does not exist.'); end
    if any((el(:,1)-el(:,2))==0), error('Both sides of element seem connected to the same node.'); end
    
    %--plot elements--
    eh = [];
    for i=1:nel
        % Start (xp) & end (xq) of elements
        xp = no(el(i,1),:);
        xq = no(el(i,2),:);

        %check element length: if not too short, leave some space for the
        %element number tag
        %{
        minlen = ;
        d = 0.1;
        len = norm(xq - xp);
        if len > minlen
            x1 = xp+1/len*(xq-xp)*(len/2-d/2);
            x2 = xp+1/len*(xq-xp)*(len/2+d/2);
            eh(end+1) = plot3(...
                [xp(1) x1(1)],[xp(2) x1(2)],[xp(3) x1(3)]);
            eh(end+1) = plot3(...
                [x2(1) xq(1)],[x2(2) xq(2)],[x2(3) xq(3)]);
                
        else
            %use index end+1 instead of i, because elements have one or two
            %line handles, depending on element length
                eh(end+1) = plot3(...
                    [xp(1) xq(1)],[xp(2) xq(2)],[xp(3) xq(3)]);
        end
%}
        
        % Plot elements
        eh(i) = plot3([xp(1) xq(1)],[xp(2) xq(2)],[xp(3) xq(3)]);
    end



    %--plot nodes--
    nh = [];
    for i=1:nno
        nh(i) = plot3(no(i,1),no(i,2),no(i,3));
    end
    
    %get largest and characteristic distance of system
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
    [normalvec, mech_dim, normalvec2] = getPlane(no); %normalvec2 only used for 1-D mechanisms, a specialty case
        
    %--plot node numbers--
    nnh = gobjects(nno,1);
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
    
    %--legend style--
    if showlegend
        set(leg1,'Marker','.','MarkerSize',16,'Color','k')
        set(leg2,'Marker','.','MarkerSize',16,'Color','g')
        set(leg3,'Marker','.','MarkerSize',16,'Color','b')
        set(leg4,'Marker','.','MarkerSize',16,'Color','r')
        set(leg5,'BackgroundColor','w','Color','b')
        set(leg6,'BackgroundColor','w','EdgeColor','k','FontWeight','bold')
    end
    
    %--visual properties--
    if nargin>2
    for i=1:nno
        if (i <= size(nprops,2) && isfield(nprops,'fix') && ~isempty(nprops(i).fix) && nprops(i).fix == true)
            set(nh(i),'Marker','.','MarkerSize',16,'Color','r')
        else
            set(nh(i),'Marker','.','MarkerSize',16,'Color','k')
        end
    end
    end
    
    if nargin>3
    for i=1:size(eprops,2)
        if isempty(eprops(i).flex)
            %rigid element
            set(eh(eprops(i).elems),'Color','k');
        else
            %flexible element
            set(eh(eprops(i).elems),'Color','b');
        end
    end
    
    el_eprops = cell2mat({eprops.elems}); %elements with associated eprops
    el_all = 1:nel; %all element numbers
    el_no_eprops = setdiff(el_all,el_eprops); %elements without eprops
    set(eh(el_no_eprops),'Color',0.75*[1 0 0]);
    end
    
    %set(eh,'Color',0.75*[1 1 1],'LineWidth',2.5)
    set(nnh,'BackgroundColor','w','Color','b')
    grid on
    
    %--axes labels and projection--
    xlabel('x');
    ylabel('y');
    zlabel('z');
    axis(getAxesLim(no,arrowlength))
    set(ah,'DataAspectRatio',[1 1 1])
    if mech_dim == 3
        set(ah,'Projection','perspective')
    end
    
    %--set view--
    try
        if mech_dim == 1
            [azi,ele] = view(normalvec);
            [azi2,ele2] = view(normalvec2);
            azis = [azi azi2];
            eles = [ele ele2];
            ii = [ele ele2]>=0;
            switch sum(ii)
                case 0
                    view(azi,ele)
                case 1
                    view(azis(ii),eles(ii))
                case 2
                    [~,mi] = min(eles);
                    view(azis(mi),eles(mi));
            end
        else
            [azi,ele] = view(normalvec);
            %do not view mechanism from below horizon, just override
            if ele<3
                %TO DO: for 2-D mechanisms, this should be different
                view(azi,15)
            end
        end
        
    catch    
        disp('Issue setting viewport.')
    end
        
%     --plot force--
%     for i=1:size(nopr,2)
%         if(isfield(nopr(i),'force') && ~isempty(nopr(i).force))
%             arrow3(no(i,:)-normVec(nopr(i).force)*arrowlength, ...
%                     no(i,:)-normVec(nopr(i).force)*clen/2,[],1.5,3)
%         end
%     end

    %---element numbers---
    %this should happen after axis properties are set because of conversions
    %between data and normalized (figure-bound) units
    %get position of labels for elements
    p(:,1:3) = (no(el(:,1),:) + no(el(:,2),:))/2;
    
    %check for coinciding coordinates
    [~,~,ic] = unique(p,'rows','stable'); %unique rows and indices
    [n,~] = histcounts(ic,1:nel); %counts of duplicates
    dupl = find(n>1); %first index of duplicate row
    
    for i = dupl
        ii = find(ic==ic(dupl)); %indices of all rows with this particular same value
        p(ii,1:3) = no(el(ii,1),:) + (no(el(ii,2),:)-no(el(ii,1),:))/3;
    end
    
    %store in guidata so accessible from other functions
    data.p = p;
    data.enh = [];
    guidata(fh,data);
    
    %plot element numbers
    %in separate function so it can also be called as a post-rotation callback
    updateElemLabel(fh); 
    
    %---end element numbers---
    
end

function updateElemLabel(fig)
    
    data = guidata(fig);
    p = data.p;
    enh = data.enh;
    
    if isempty(enh)
        %plot
        enh = gobjects(size(p,1),1);
        for i=1:size(p,1)
            enh(i) = text(p(i,1),p(i,2),p(i,3),num2str(i),'Units','data','visible','off');
            set(enh(i),'Units','normalized')
            set(enh(i),'Visible','on')
            set(enh(i),'BackgroundColor','w','EdgeColor','k','FontWeight','bold');
        end
        %store handles
        data.enh = enh;
        guidata(fig,data);
    else
        %just update position
        for i=1:size(p,1)
            set(enh(i),'Units','data')
            set(enh(i),'Position',p(i,1:3))
            set(enh(i),'Units','normalized')
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
    no_c = no - repmat(mean(no,1),size(no,1), 1);
    
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

function rotationcallback(~,~)
    %normally, inputs are obj and evd
    
    updateElemLabel(gcbo);
    %gcbo contains handle of figure that caused callback
    
end
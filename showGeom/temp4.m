function temp4

    close
    clear
    clc

    fh = figure;
    hold on

    ph = plot3([0 1],[0 1],[0 1]);

    set(gcf,'Position',[19   285   560   420])
    set(gca,'view',[9.2000   23.6000])
    grid on

    k = rotate3d;
    k.ActionPostCallback = @mypostcallback;
    k.Enable = 'on';

    p = [0.2 0.2 0.2;
        0.3 0.3 0.3];
    
	data.p = p;
    data.elh = [];
    guidata(fh,data);
    
    scatter3(p(:,1),p(:,2),p(:,3));
    
    updateElemLabel(fh);
    
end

function updateElemLabel(fig)
    
    data = guidata(fig);
    p = data.p;
    elh = data.elh;
    
    if ~isempty(elh)
        delete(elh)
    end
    
    elh = [];
    for i=1:size(p,1)
        temph = text(p(i,1),p(i,2),p(i,3),num2str(i),'Units','data','visible','off');
        temph.Units = 'normalized';
        np = temph.Position;
        nx = np(1);
        ny = np(2);
        delete(temph);
        elh(i) = text(nx,ny,num2str(i),'Units','normalized');
        set(elh(i),'BackgroundColor','white');
    end
    data.elh = elh;
    guidata(fig,data);
    
end

function mypostcallback(obj,evd)
    
    updateElemLabel(gcbo); 
    %gcbo contains handle of figure that caused callback   
    
end
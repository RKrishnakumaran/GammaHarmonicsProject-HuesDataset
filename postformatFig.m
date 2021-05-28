% figno = 1;
% f = openfig(['' num2str(floor(figno))]);
% destfolder = ['' num2str(floor(figno))]
% destfilename = ['' num2str(floor(figno))]
%%
f;
for i = numel(f.Children):-1:1
    x = f.Children(i); 
    if strcmp( x.Type, 'legend')
        set(x, 'FontSize', 12); set(x,'Linewidth',1);
    elseif strcmp( x.Type, 'axes') 
        box off; set(x,'TickDir','out'); set(x,'TickLength',[0.02, 0.05]); set(x,'Linewidth',1); set(x,'FontSize',16);
    elseif strcmp(x.Type, 'polaraxes')
        box off; set(x,'TickDir','out'); set(x,'TickLength',get(x,'TickLength')); set(x,'Linewidth',1); set(x,'FontSize',16);
    end
for j = 1:numel(x.Children)
    y = x.Children(j);
    if strcmp( y.Type, 'text')
        set(y, 'FontSize', 14);
    elseif strcmp( y.Type, 'line')
        
        if strcmp( y.Marker, 'none')
            set(y, 'LineWidth', 1.5);
        else
            disp('higha')
            set(y, 'MarkerSize', 18);
        end
    elseif strcmp( y.Type, 'scatter')
        set(y, 'SizeData', 100);
    end
end
end

%% Correcting marker color in fig 1 last subplot
f;
for i = 1:numel(f.Children)
    x = f.Children(i);
    allchildx = allchild(x);
    for j = 1:numel(allchildx)
        y = allchildx(j);
        if strcmp( y.Type, 'line')
            if ~strcmp( y.Marker , 'none')
                set(y, 'color', 'r');
                set(y,'Marker','x');
                set(y,'linewidth', 1);
            end
        end
    end
end
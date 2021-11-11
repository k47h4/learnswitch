% This script takes the data from the local directory and makes
% panels for Fig 6B,C and D


do_ResponseCurves = 1;
do_BarPlots = 1;
do_Matrix = 1;
do_Square = 1;

if do_ResponseCurves % reponse curves
    D1 = load(['Fig6B_responses.mat']);
    
    D1.time
    D1.time_unit
    D1.means
    
    figure;
    plot(D1.time, D1.means.PYR,'k');
    hold on;plot(D1.time, D1.means.PV,'Color',[255 165 0]./255);
    hold on;plot(D1.time, D1.means.SST,'b');
    hold on;plot(D1.time, D1.means.VIP,'g');
    legend({'PYR','PV','SOM','VIP'});
end


%-------------------------------------------------------------------------
if do_BarPlots % bar plots
    D2 = load(['Fig6D']);
    
    D2.peaks
    D2.xlabels
    
    figure;bar([D2.peaks.PYR,D2.peaks.PV,D2.peaks.SOM,D2.peaks.VIP]);
end

%-------------------------------------------------------------------------

if do_Matrix % All manipulations matrix
    fn = {
        'rsc_multiplicative.mat' 
        'SI_multiplicative.mat'     
        'rsc_additive.mat' 
        'SI_additive.mat'}; 
    
    for i=1:length(fn) 
        
        D = load([fn{i}]);
        
        % re-order
        reorderedY =  [1, 4, 2, 3, 5:15];
        if i==1||i==3
            reorderedX =  [3,10,7,1,9,2,5,6,4,8];
        else
            reorderedX =  [2, 1, 3, 4];
        end
        
        D.changes = D.changes(reorderedY,reorderedX);
        D.xlabels2 = D.ylabels(reorderedX,:);
        D.ylabels = D.xlabels(reorderedY,:);
        D.xlabels = D.xlabels2;
        
        %remove _
        for s = 1:size(D.xlabels,1)
            D.xlabels(s,:)=strrep(D.xlabels(s,:),'_','-');
        end
        for s = 1:size(D.ylabels,1)
            D.ylabels(s,:)=strrep(D.ylabels(s,:),'_','-');
        end
        
        figure;imagesc(D.changes);
        colorbar;colormap(linspecer)
        c=colormap;cL = round(size(c,1)/2);
        set(gca,'YTick',[1:size(D.changes,1)],'YTickLabel',D.ylabels);
        set(gca,'XTick',[1:size(D.changes,2)],'XTickLabel',D.xlabels);
        xtickangle(90)
        
        % grids
        [rows, columns] = size(D.changes);
        hold on;
        for row = 1 : rows+1
            line([-.5, columns+.5], [row-.5, row-.5], 'Color', 'w','linewidth',1);
        end
        for col = 1 : columns+1
            line([col-.5, col-.5], [-.5, rows+.5], 'Color', 'w','linewidth',1);
        end
        
        ext='-dpdf';fext='.pdf';
        if i==1
            set(gca,'CLim',[-0.5,0.5]);
            set(gcf,'position',[140  250  490  420])
            filename = sprintf('%s%s_ROIsImage%s',rdir,'matrixNCmult',fext);
            title('matrixNCmult')
        elseif i==2
            colormap(c(cL:end,:))
            set(gca,'CLim',[0,240]);
            set(gcf,'position',[635  250  275  420])
            filename = sprintf('%s%s_ROIsImage%s',rdir,'matrixSELmult',fext);
            title('matrixSELmult')
        elseif i==3
            set(gca,'CLim',[-0.5,0.5]);
            set(gcf,'position',[140  250  490  420])
            filename = sprintf('%s%s_ROIsImage%s',rdir,'matrixNCadd',fext);
            title('matrixNCadd')
        elseif i==4
            colormap(c(cL:end,:))
            set(gca,'CLim',[0,100]);
            set(gcf,'position',[635  250  275  420])
            filename = sprintf('%s%s_ROIsImage%s',rdir,'matrixSELadd',fext);
            title('matrixSELadd')
        end
        %hf=gcf;eval(sprintf('print -f%d %s ''%s''',hf.Number,ext,filename))
    end
end


%----------------------------------------------------------

if do_Square % square plot
    fn = {
    'Fig6E.mat'   
    };
    ext='-dpdf';fext='.pdf';
    for i=1:length(fn) % i=1 i=2
        uM2 = [];
    
        D = load([fn{i}]);
    
        xIgnore = find (D.SOM_modulation == 1);
        xAttend = find (D.SOM_modulation == 2.2);
        l=fieldnames(D.rsc);
    
        for j=1:length(l)
            y=D.rsc.(l{j});
            uM2 = [uM2;y(xIgnore) y(xAttend)];
        end
        square_plot (uM2);
        t = strrep(fn{i},'_','-');t=t(1:end-4);
        title(t)
        filename = sprintf('%s%s%s',rdir,t,fext);
        %hf=gcf;eval(sprintf('print -f%d %s ''%s''',hf.Number,ext,filename))
    end

end

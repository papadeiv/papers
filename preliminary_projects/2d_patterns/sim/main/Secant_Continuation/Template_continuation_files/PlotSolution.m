function plotHandle = PlotSolution(u,p,parentHandle)

   if isempty(parentHandle)
     scrsz = get(0,'ScreenSize');
     plotHandle = figure('Position',[2/4*scrsz(3) scrsz(4)/2 scrsz(3)/4 scrsz(4)/4]);
     parentHandle = plotHandle;
   else
     plotHandle = parentHandle;
   end 

   figure(parentHandle);
   
     % Rename parameters

  % plot
drawnow;


   print -dtiff state.tiff

end

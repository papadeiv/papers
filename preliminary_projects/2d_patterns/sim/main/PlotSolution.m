function plotHandle = PlotSolution(uu,p,parentHandle,mesh_params)

   if isempty(parentHandle)
     scrsz = get(0,'ScreenSize');
     plotHandle = figure('Position',[2/4*scrsz(3) scrsz(4)/2 scrsz(3)/4 scrsz(4)/4]);
     parentHandle = plotHandle;
   else
     plotHandle = parentHandle;
   end 

   figure(parentHandle);
   
   % Extract the variables
   N = mesh_params.N; 
   r = mesh_params.r;
   u_re = uu(1:N);
   u_im = uu(N+1:2*N);
  
   % Produce the plot
   plot(r,u_re,'b',r,u_im,'r');  drawnow;
   %print -dtiff state.tiff

end

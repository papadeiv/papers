function plotHandle = PlotSolution_2d(uu,mesh_params)
   
   % Extract the variables
   N = mesh_params.N; 
   r = mesh_params.r;
   u_re = uu(1:N);
   u_im = uu(N+1:2*N);

   % Angular coordinate (180°)
    nTheta = 500;
    theta = linspace(0, 2*pi, nTheta); % 0 to 180 degrees

    % Polar grid    
    [R, Theta] = meshgrid(r, theta);

    % Replicate radial profile along angular direction
    F = repmat(u_re', nTheta, 1);

    % Convert to Cartesian coordinates
    X = R .* cos(Theta);
    Y = R .* sin(Theta);

    % Filled contour plot
    figure;
    contourf(X, Y, F, 20, 'LineColor', 'none'); % 20 contour levels
    axis equal tight
    colormap jet
    colorbar
    xlabel('x'); ylabel('y');
    title('Radially symmetric contour plot');
    drawnow;

end

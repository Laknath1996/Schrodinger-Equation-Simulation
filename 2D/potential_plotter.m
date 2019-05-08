function potential_plotter(x,y,V)
    figure('Name','Potential','position',[100 100 1200 400])
   
    subplot(1,2,1)
    hold on;grid on;
    surf(x*1e10,y*1e01,V)

    colormap(jet)
    colorbar
    view(30,30)
    %shading flat

    xlabel('x (angstroms)')
    ylabel('y (angstroms)')
    zlabel('Energy (eV)')


    subplot(1,2,2)
    hold on;grid on;
    pcolor(x*1e10,y*1e10,V)

    colormap(jet)
    colorbar
    %shading flat

    xlabel('x (angstroms)')
    ylabel('y (angstroms)')
    zlabel('Energy (eV)')
end

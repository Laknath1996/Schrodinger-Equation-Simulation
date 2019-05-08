function wavefunction_plotter(x,y,V,num_sol,psi1,E)

    c=0;
    ii=0;
    for i=1:num_sol
        if i>45
          break
        end
        if i==1 || i==16 || i==31
          figure('Name','FEM method','position',[100 100 1600 900])
          c=c+1;
          ii=0;
        end
        ii=ii+1;

        subplot(3,5,ii,'fontsize',10)
        hold on

        pcolor(x*1e10,y*1e10,(psi1(:,:,i)) )  
        contour(x*1e10,y*1e10,V,1,'linewidth',3,'linecolor','w')

        xlabel('x (angstroms)')
        ylabel('y (angstroms)')
        title(strcat('E',num2str(num_sol+1-i),'=',num2str(E(i,1)*1000,'%.1f'),'eV'))
        %axis equal
        shading flat
        colormap(jet)
    end

end

function energy_plotter(ee,z,V)
    figure;
    plot(z*10^10,V,'r-','LineWidth',1.5);ylabel('V_o (eV)','fontweight','bold');xlabel('z (Angstroms)','fontweight','bold');
    hold on;
    for i = 1:length(ee)
        line([min(z*10^10) max(z*10^10)],[ee(i) ee(i)],'Color','b','LineWidth',1.5);
        text(0,ee(i)-0.25,sprintf('%.3f eV',ee(i)))
    end
    title('Allowed Energies');
    legend('V(z)','Allowed E');
end

        
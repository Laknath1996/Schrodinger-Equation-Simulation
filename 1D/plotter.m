function plotter(ee,z,V,dz,N)
    % % plot the wavefunctions
    for i = 1:length(ee)
        psi = wav_func(ee(i),V,dz,N);
        figure;
        subplot(211);
        yyaxis left
        plot(z*10^10,psi,'LineWidth',1.1);xlabel('z (Angstroms)','fontweight','bold');ylabel('\psi','fontweight','bold','fontsize',16);hold on;
        yyaxis right
        plot(z*10^10,V,'r-','LineWidth',1.5);ylabel('V_o (eV)','fontweight','bold');
        title(sprintf('Wavefunction for n = %d',i));
        subplot(212);
        yyaxis left
        plot(z*10^10,psi.^2,'k','LineWidth',1.1);xlabel('z (Angstroms)','fontweight','bold');ylabel('\psi^2','fontweight','bold','fontsize',16);hold on;
        yyaxis right
        plot(z*10^10,V,'r-','LineWidth',1.5);ylabel('V_o (eV)','fontweight','bold');
        title(sprintf('P(z) for n = %d',i));

    end
end
calk_initpy( ...
    'C:\Users\yau17reu\anaconda\Anaconda3\envs\calkenv\pythonw.exe')

[Macid,pH,Tk,Msamp,Cacid,S,XT,KX] = calk_Dickson1981;


figure(1); clf

printsetup(gcf, [10 7])

scatter(Macid*1e3,pH,10,'filled')
xlabel('Acid mass / g')
ylabel('Free scale pH')
grid on

setaxes(gca,10)
set(gca, 'xticklabel',num2str(get(gca,'xtick')','%.1f'))

ylim([2 9])

% print('-r300','docs/figures/Macid_pH_D81','-dpng')
function t = calk_ptl(datpath,datfile,Vsamp,Cacid,CT,S,PT,SiT, ...
    burette_cx,Tk_force,printpath)

% Get initial Gran guesses
clear t
[t.Macid,t.EMF,t.Tk,Msamp, t.F1g,t.Lg,t.EMF0g, ~,EMF0g,t.pHg,t.L] ...
    = calk_guessGran([datpath datfile],Vsamp,Cacid,S);
t = struct2table(t);

t.mu = Msamp ./ (Msamp + t.Macid);

% Solve for TA fully
[AT,EMF0,AT_RMS,AT_Npts] = calk_VINDTA([datpath datfile],Vsamp,Cacid, ...
    S,CT,PT,SiT,burette_cx,Tk_force);

[t.H,t.pH] = calk_EMF2H(t.EMF,EMF0,t.Tk);

[t.simAT,components] = calk_simAT(t.Macid,t.Tk,t.H,Msamp,S,CT,PT,SiT);
cvars = fieldnames(components);
t = [t struct2table(components)];

t.estAT = (t.simAT + t.Macid*Cacid ./ (t.Macid + Msamp)) ./ t.mu;

% ----------------------------------------------------- Plot the lot! -----

fvars = {'EMF' 'estAT' 'F1g' [] 'EMF0g'};

fclr_guess = [0.96 0.86 0.04];
fclr_final = [0.21 0.46 1];
fclr_both = [0.27 0.8 0.54];
mksz = 20;

ffsz = 8;

fxlim = minmax(1e3*t.Macid');

freg = regstats(1e-6*t.F1g(t.Lg),1e3*t.Macid(t.Lg),'linear','beta');
fxint = -freg.beta(1)/freg.beta(2);

ftxt_mu     = char(956);
ftxt_endash = char(8211);

figure(1); clf

set(gcf, 'color','w', 'paperunits','centimeters', 'units','centimeters')
set(gcf, 'papersize',[18 18])
set(gcf, 'paperposition',[0 0 get(gcf,'papersize')])

for V = 1:numel(fvars)
if ~isempty(fvars{V})
    
subplot(3,2,V); hold on

xlim(fxlim)

switch fvars{V}

    case 'EMF'

        ylim(minmax(t.(fvars{V})')+30*[-1 1])
        set(gca, 'ytick',0:100:800)
        
        VLg = t.Lg & ~t.L;
        VL  = t.L & ~t.Lg;
        VLb = t.Lg & t.L;
        VLn = ~t.Lg & ~t.L;
        
        nl = plot(1e3*AT*Msamp/Cacid*[1 1],get(gca,'ylim'), ...
            'color',fclr_final, 'linewidth',1); calk_nolegend(nl)
        
        nl = plot(fxint*[1 1],get(gca,'ylim'), ...
            'color',fclr_guess, 'linestyle','--', 'linewidth',1); 
        calk_nolegend(nl)
        
        scatter(1e3*t.Macid(VLg),t.(fvars{V})(VLg),mksz,fclr_guess, ...
            'filled', 'markeredgecolor','k')
        scatter(1e3*t.Macid(VL) ,t.(fvars{V})(VL) ,mksz,fclr_final, ...
            'filled', 'markeredgecolor','k')
        scatter(1e3*t.Macid(VLb),t.(fvars{V})(VLb),mksz,fclr_both, ...
            'filled', 'markeredgecolor','k')
        scatter(1e3*t.Macid(VLn),t.(fvars{V})(VLn),mksz,'k')
        
        ylabel('EMF / mV')
        
        legend('First guess','Final result','Both', 'location','nw')
        
        text(0,1.1,['(a) Final EMF° = ' num2str(EMF0,'%.2f') ...
            ' mV'], 'fontname','arial', ...
            'fontsize',ffsz, 'color','k', 'units','normalized')

    case 'F1g'
        
        ylim([0 max(1e-6*t.(fvars{V}))*1.1])
        set(gca, 'ytick',0:4:50)

        ylabel(['\itF\rm_1 \times 10^{' ftxt_endash '6}'])
        
        plot(fxint*[1 1],get(gca,'ylim'), ...
            'color',fclr_guess, 'linestyle','--', 'linewidth',1)
        
        plot(fxlim,x2fx(fxlim','linear')*freg.beta, 'color',fclr_guess, ...
            'linewidth',1)
        
        scatter(1e3*t.Macid(t.Lg),1e-6*t.(fvars{V})(t.Lg),mksz, ...
            fclr_guess,'filled', 'markeredgecolor','k')
        scatter(1e3*t.Macid(~t.Lg),1e-6*t.(fvars{V})(~t.Lg),mksz,'k')
        
        text(0,1.1,['(b) First guess \itA\rm_T = ' ...
            num2str(1e3*fxint*Cacid/Msamp,'%.1f') ...
            ' ' ftxt_mu 'mol\cdotkg^{' ftxt_endash '1}'], ...
            'fontname','arial', 'fontsize',ffsz, ...
            'color','k', 'units','normalized')
        
        
    case 'EMF0g'
        
        ylim(minmax(t.(fvars{V})') .* [0.999 1.001])
        set(gca, 'ytick',600:700)
        
        ylabel('First guess EMF° / mV')
        
        plot(fxint*[1 1],get(gca,'ylim'), ...
            'color',fclr_guess, 'linestyle','--', 'linewidth',1)
        
        plot(fxlim,EMF0g*[1 1], 'color',fclr_guess, 'linewidth',1)

        scatter(1e3*t.Macid(t.Lg),t.(fvars{V})(t.Lg),mksz,fclr_guess, ...
            'filled', 'markeredgecolor','k')
        scatter(1e3*t.Macid(~t.Lg),t.(fvars{V})(~t.Lg),mksz,'k')

        text(0,1.1,['(c) First guess EMF° = ' num2str(EMF0g,'%.1f') ...
            ' mV'], 'fontname','arial', 'fontsize',ffsz, ...
            'color','k', 'units','normalized')
        
        
    case 'estAT'
        
        plot(fxlim,AT*[1 1]*1e6, 'color',fclr_final, 'linewidth',1)
        
        scatter(1e3*t.Macid(t.L),t.(fvars{V})(t.L)*1e6,mksz,fclr_final, ...
            'filled', 'markeredgecolor','k')
        scatter(1e3*t.Macid(~t.L),t.(fvars{V})(~t.L)*1e6,mksz,'k')
        
        ylim(minmax(t.(fvars{V})') .* [0.995 1.005] * 1e6)
        
        ylabel(['\itA\rm_T from pH / ' ftxt_mu 'mol\cdotkg^{' ...
            ftxt_endash '1}'])
        
        text(0,1.1,['(d) Final \itA\rm_T = ' num2str(1e6*AT,'%.1f') ...
            ' ± ' num2str(1e6*AT_RMS,'%.1f') ' ' ftxt_mu ...
            'mol\cdotkg^{' ftxt_endash '1} (\itn\rm = ' ...
            num2str(AT_Npts) ')'], ...
            'fontname','arial', 'fontsize',ffsz, ...
            'color','k', 'units','normalized')
        

end %switch

xlabel('Acid mass / g')

set(gca, 'box','on', 'tickdir','out', 'xcolor','k', 'ycolor','k', ...
    'fontname','arial', 'fontsize',ffsz, 'xtick',0:10, ...
    'labelfontsizemultiplier',1)

end %if
end %for V

subplot(3,2,[4 6]); hold on

cvars_names = { ...
    ['+[HCO_3^' ftxt_endash ']'] ...
    ['+2[CO_3^{2' ftxt_endash '}]'] ...
    ['+[B(OH)_4^' ftxt_endash ']'] ...
    ['+[OH^' ftxt_endash ']'] ...
    [ftxt_endash '[HSO_4^' ftxt_endash ']'] ... ...
    [ftxt_endash '[HF]'] ...
    [ftxt_endash '[H_3PO_4]'] ...
    ['+[HPO_4^{2' ftxt_endash '}]'] ...
    ['+2[PO_4^{3' ftxt_endash '}]'] ...
    ['+[SiO(OH)_3^' ftxt_endash ']'] ...
    [ftxt_endash '[H^+]']};
cvars_clrs = { ...
    0.4*[1 1 1] ...
    0.4*[1 1 1] ...
    [0.19 0.31 0.97] ...
    [1 0.05 0.05] ...
    [0.95 0.8 0.14] ...
    [0.12 0.74 0.12] ...
    [1 0.5 0] ...
    [1 0.5 0] ...
    [1 0.5 0] ...
    [0.94 0.56 0.63] ...
    [1 0.05 0.05]};

fcvars = [cvars; 'H'];

clr_AT = 'k';

    xlim(fxlim)
    ylim([1e-16 1e-2])
    set(gca, 'ytick',10.^(-16:2:4))
    
    fpxL = t.pHg >=3 & t.pHg <= 4;
    fpx = [min(t.Macid(fpxL)) max(t.Macid(fpxL))]*1e3;
    spyl = get(gca,'ylim');
    patch(fpx([1 2 2 1 1]),spyl([1 1 2 2 1]),fclr_final, ...
        'edgecolor','none', 'facealpha',0.2)

    plot(1e3*t.Macid,t.estAT, 'color',clr_AT, 'linewidth',1, ...
        'marker','o', 'markersize',2)
    text(max(1e3*t.Macid)*1.02,t.estAT(end),'\itA\rm_T', ...
        'fontname','arial', 'fontsize',ffsz, 'color',clr_AT)

    for V = 1:numel(fcvars)

        plot(1e3*t.Macid,abs(t.(fcvars{V})), ...
            'color',cvars_clrs{V}, 'marker','o', 'markersize',2)
        text(max(1e3*t.Macid)*1.02,abs(t.(fcvars{V})(end)), ...
            cvars_names{V}, 'fontname','arial', 'fontsize',ffsz, ...
            'color',cvars_clrs{V})

    end %for V

%     legend(['estAT'; cvars], 'location','eastoutside')

    set(gca, 'box','on', 'tickdir','out', 'xcolor','k', 'ycolor','k', ...
        'fontname','arial', 'fontsize',ffsz, 'xtick',0:10, ...
        'yscale','log', 'labelfontsizemultiplier',1)
    
    xlabel('Acid mass / g')
    ylabel(['Concentration from pH / mol\cdotkg^{' ftxt_endash '1}'])
    
    text(0,1.04,'(e) \itA\rm_T components', 'fontname','arial', ...
        'fontsize',ffsz, 'color','k', 'units','normalized')

    
annotation('textbox',[0.25 0.9 0.5 0.1], 'string',datfile, ...
    'fontname','arial', 'fontsize',ffsz*1.1, 'edgecolor','none', ...
    'fontweight','bold', 'horizontalalignment','center')

datnodat = datfile(1:end-4);
if ~isempty(printpath)
    print('-r300',[printpath 'calk_ptl_' datnodat],'-dpng')
end %if
    
end %function calk_ptl

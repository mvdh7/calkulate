function [t,AT,EMF0,AT_RMS,AT_Npts] = calk_ptl(datpath,datfile,...
    Vsamp,Cacid,CT,S,PT,SiT,burette_cx,Tk_force,printpath)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% calk_ptl ("plot the lot") draws - hopefully - eveything you would ever
%  want to know about a single total alkalinity titration.
% 
% Part of https://github.com/mvdh7/calkulate
% Written by Matthew P. Humphreys, last updated 2018-01-15
% 
% ~~~~~ INPUTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%     datpath: string, 'path/to/datfiles/'
%     datfile: string, 'dat file name.dat'
%       Vsamp: sample volume in ml
%       Cacid: acid concentration in mol/kg (e.g. from calk_VINDTA_CRM)
%          CT: dissolved inorganic carbon in mol/kg-sw
%           S: practical salinity
%          PT: total phosphate in mol/kg-sw
%         SiT: total silicate in mol/kg-sw
%  burette_cx: burette correction factor (1 if unknown)
%    Tk_force: override titration file temperature in K ([] if not needed)
%   printpath: path to print figure to ([] for no figure)
%
% ~~~~~ OUTPUTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%           t: table of titration data
%          AT: final total alkalinity in mol/kg-sw
%        EMF0: final EMF0 in mV
%      AT_RMS: root-mean-square of AT estimates in mol/kg-sw
%     AT_Npts: number of data points used for final AT estimate
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


% ====================================================== Process data =====

% Get initial Gran guesses
clear t
[t.Macid,t.EMF,t.Tk,Msamp, t.F1g,t.Lg,t.EMF0g, ~,EMF0g,t.pHg,t.L] ...
    = calk_guessGran([datpath datfile],Vsamp,Cacid,S);
t.mu = Msamp ./ (Msamp + t.Macid);
t = struct2table(t);

% Solve for TA fully and get components
[AT,EMF0,AT_RMS,AT_Npts] = calk_VINDTA([datpath datfile],Vsamp,Cacid, ...
    S,CT,PT,SiT,burette_cx,Tk_force);

[t.H,t.pH] = calk_EMF2H(t.EMF,EMF0,t.Tk);

[t.simAT,components] = calk_simAT(t.Macid,t.Tk,t.H,Msamp,S,CT,PT,SiT);
cvars = fieldnames(components);
t = [t struct2table(components)];

t.estAT = (t.simAT + t.Macid*Cacid ./ (t.Macid + Msamp)) ./ t.mu;

% Perform first guess alkalinity regression
freg = regstats(1e-6*t.F1g(t.Lg),1e3*t.Macid(t.Lg),'linear','beta');
fxint = -freg.beta(1)/freg.beta(2);


% ================================== Define settings and begin figure =====

% Colours for first guess/final result highlights
fclr_guess = [0.96 0.86 0.04];
fclr_final = [0.21 0.46 1];
fclr_both = [0.27 0.8 0.54];

% Marker sizes
mksz = 20;

% Font size
ffsz = 8;

% x-axis limits
fxlim = minmax(1e3*t.Macid');

% Define special characters
ftxt_mu     = char(956);
ftxt_endash = char(8211);

% Begin figure
figure(1); clf

% Figure print settings
set(gcf, 'color','w', 'paperunits','centimeters', 'units','centimeters')
set(gcf, 'papersize',[18 18])
set(gcf, 'paperposition',[0 0 get(gcf,'papersize')])


% ==================================================== Begin subplots =====


% -------------------------------------------------- (a) Measured EMF -----

subplot(3,2,1); hold on

    % Limits
    xlim(fxlim)
    ylim(minmax(t.EMF')+30*[-1 1])
    
    % Logicals
    VLg = t.Lg & ~t.L;
    VL  = t.L & ~t.Lg;
    VLb = t.Lg & t.L;
    VLn = ~t.Lg & ~t.L;

    % Final alkalinity
    nl = plot(1e3*AT*Msamp/Cacid*[1 1],get(gca,'ylim'), ...
        'color',fclr_final, 'linewidth',1); calk_nolegend(nl)

    % First guess alkalinity
    nl = plot(fxint*[1 1],get(gca,'ylim'), ...
        'color',fclr_guess, 'linestyle','--', 'linewidth',1); 
    calk_nolegend(nl)

    % EMF measurements
    scatter(1e3*t.Macid(VLg),t.EMF(VLg),mksz,fclr_guess, ...
        'filled', 'markeredgecolor','k')
    scatter(1e3*t.Macid(VLb),t.EMF(VLb),mksz,fclr_both, ...
        'filled', 'markeredgecolor','k')
    scatter(1e3*t.Macid(VL) ,t.EMF(VL) ,mksz,fclr_final, ...
        'filled', 'markeredgecolor','k')
    scatter(1e3*t.Macid(VLn),t.EMF(VLn),mksz,'k')

    % Labels
    xlabel('Acid mass / g')
    ylabel('EMF / mV')

    legend('First guess','Both','Final result', 'location','nw')

    text(0,1.1,['(a) Final EMF° = ' num2str(EMF0,'%.2f') ...
        ' mV'], 'fontname','arial', ...
        'fontsize',ffsz, 'color','k', 'units','normalized')
    
    % Axis settings
    set(gca, 'box','on', 'tickdir','out', 'xcolor','k', 'ycolor','k', ...
        'fontname','arial', 'fontsize',ffsz, 'xtick',0:10, ...
        'labelfontsizemultiplier',1, 'ytick',0:100:800)

    
% ---------------------------------------- (b) First guess alkalinity -----

subplot(3,2,3); hold on

    % Limits
    xlim(fxlim)
    ylim([0 max(1e-6*t.F1g)*1.1])
 
    % First guess alkalinity
    plot(fxint*[1 1],get(gca,'ylim'), ...
        'color',fclr_guess, 'linestyle','--', 'linewidth',1)

    % First guess alkalinity regression line
    plot(fxlim,x2fx(fxlim','linear')*freg.beta, 'color',fclr_guess, ...
        'linewidth',1)

    % F1 function points
    scatter(1e3*t.Macid(t.Lg),1e-6*t.F1g(t.Lg),mksz, ...
        fclr_guess,'filled', 'markeredgecolor','k')
    scatter(1e3*t.Macid(~t.Lg),1e-6*t.F1g(~t.Lg),mksz,'k')
    
    % Labels
    xlabel('Acid mass / g')
    ylabel(['\itF\rm_1 \times 10^{' ftxt_endash '6}'])
    
    text(0,1.1,['(b) First guess \itA\rm_T = ' ...
        num2str(1e3*fxint*Cacid/Msamp,'%.1f') ...
        ' ' ftxt_mu 'mol\cdotkg^{' ftxt_endash '1}'], ...
        'fontname','arial', 'fontsize',ffsz, ...
        'color','k', 'units','normalized')

    % Axis settings
    set(gca, 'box','on', 'tickdir','out', 'xcolor','k', 'ycolor','k', ...
        'fontname','arial', 'fontsize',ffsz, 'xtick',0:10, ...
        'labelfontsizemultiplier',1, 'ytick',0:4:50)
        
        
% ---------------------------------------------- (c) First guess EMF0 -----

subplot(3,2,5); hold on
        
    % Limits
    xlim(fxlim)
    ylim(minmax(t.EMF0g') .* [0.999 1.001])

    % First guess alkalinity
    plot(fxint*[1 1],get(gca,'ylim'), ...
        'color',fclr_guess, 'linestyle','--', 'linewidth',1)

    % First guess EMF0
    plot(fxlim,EMF0g*[1 1], 'color',fclr_guess, 'linewidth',1)

    % First guess EMF0 points
    scatter(1e3*t.Macid(t.Lg),t.EMF0g(t.Lg),mksz,fclr_guess, ...
        'filled', 'markeredgecolor','k')
    scatter(1e3*t.Macid(~t.Lg),t.EMF0g(~t.Lg),mksz,'k')

    % Labels
    text(0,1.1,['(c) First guess EMF° = ' num2str(EMF0g,'%.1f') ...
        ' mV'], 'fontname','arial', 'fontsize',ffsz, ...
        'color','k', 'units','normalized')
    
    xlabel('Acid mass / g')
    ylabel('First guess EMF° / mV')

    % Axis settings
    set(gca, 'box','on', 'tickdir','out', 'xcolor','k', 'ycolor','k', ...
        'fontname','arial', 'fontsize',ffsz, 'xtick',0:10, ...
        'labelfontsizemultiplier',1, 'ytick',600:700)
        
        
% ---------------------------------------- (d) Final total alkalinity -----

subplot(3,2,2); hold on

    % Limits
    xlim(fxlim)
    ylim(minmax(t.estAT') .* [0.995 1.005] * 1e6)
        
    % Final alkalinity
    plot(fxlim,AT*[1 1]*1e6, 'color',fclr_final, 'linewidth',1)

    % Alkalinity points
    scatter(1e3*t.Macid(t.L),t.estAT(t.L)*1e6,mksz,fclr_final, ...
        'filled', 'markeredgecolor','k')
    scatter(1e3*t.Macid(~t.L),t.estAT(~t.L)*1e6,mksz,'k')

    % Labels
    xlabel('Acid mass / g')
    ylabel(['\itA\rm_T from pH / ' ftxt_mu 'mol\cdotkg^{' ...
        ftxt_endash '1}'])

    text(0,1.1,['(d) Final \itA\rm_T = ' num2str(1e6*AT,'%.1f') ...
        ' ± ' num2str(1e6*AT_RMS,'%.1f') ' ' ftxt_mu ...
        'mol\cdotkg^{' ftxt_endash '1} (\itn\rm = ' ...
        num2str(AT_Npts) ')'], ...
        'fontname','arial', 'fontsize',ffsz, ...
        'color','k', 'units','normalized')
        
    % Axis settings
    set(gca, 'box','on', 'tickdir','out', 'xcolor','k', 'ycolor','k', ...
        'fontname','arial', 'fontsize',ffsz, 'xtick',0:10, ...
        'labelfontsizemultiplier',1, 'ytick',0:10:1e4)


% ----------------------------------------- (e) Alkalinity components -----
   
% Define labels and colours for the components
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

% Add hydrogen ion concentration to the list
fcvars = [cvars; 'H'];

% Define colour for total alkalinity
clr_AT = 'k';

% Begin subplot
subplot(3,2,[4 6]); hold on

    % Limits
    xlim(fxlim)
    ylim([1e-16 1e-2])
    
    % Patch to show region of data that are used
    fpxL = t.pHg >=3 & t.pHg <= 4;
    fpx = [min(t.Macid(fpxL)) max(t.Macid(fpxL))]*1e3;
    spyl = get(gca,'ylim');
    patch(fpx([1 2 2 1 1]),spyl([1 1 2 2 1]),fclr_final, ...
        'edgecolor','none', 'facealpha',0.2)

    % Total alkalinity plot
    plot(1e3*t.Macid,t.estAT, 'color',clr_AT, 'linewidth',1, ...
        'marker','o', 'markersize',2)
    text(max(1e3*t.Macid)*1.02,t.estAT(end),'\itA\rm_T', ...
        'fontname','arial', 'fontsize',ffsz, 'color',clr_AT)

    % Components of alkalinity plots
    for V = 1:numel(fcvars)
        
        VL = t.(fcvars{V}) ~= 0;
        
        if any(VL)
            
            plot(1e3*t.Macid,abs(t.(fcvars{V})), ...
                'color',cvars_clrs{V}, 'marker','o', 'markersize',2)
            text(max(1e3*t.Macid)*1.02,abs(t.(fcvars{V})(end)), ...
                cvars_names{V}, 'fontname','arial', 'fontsize',ffsz, ...
                'color',cvars_clrs{V})
            
        end %if

    end %for V

    % Axis settings
    set(gca, 'box','on', 'tickdir','out', 'xcolor','k', 'ycolor','k', ...
        'fontname','arial', 'fontsize',ffsz, 'xtick',0:10, ...
        'yscale','log', 'labelfontsizemultiplier',1, 'ytick',10.^(-16:4))
    set(gca, 'yticklabel',num2str(-log10(get(gca,'ytick'))','%.0f'))
    
    % Labels
    xlabel('Acid mass / g')
    ylabel([ftxt_endash ...
        'log_{10} (Concentration from pH / mol\cdotkg^{' ...
        ftxt_endash '1})'])
    
    text(0,1.04,'(e) \itA\rm_T components', 'fontname','arial', ...
        'fontsize',ffsz, 'color','k', 'units','normalized')

    
% ================================================ Figure annotations =====

annotation('textbox',[0.25 0.9 0.5 0.1], 'string',datfile, ...
    'fontname','arial', 'fontsize',ffsz*1.1, 'edgecolor','none', ...
    'fontweight','bold', 'horizontalalignment','center')

bottxt = ['Sample volume: ' num2str(Vsamp) ...
    ' ml :: Acid concentration: ' num2str(Cacid) ' mol/kg'];

annotation('textbox',[0 0 1 0.04], 'string',bottxt, ...
    'fontname','arial', 'fontsize',ffsz, 'edgecolor','none', ...
    'fontweight','bold', 'horizontalalignment','center')


% ====================================================== Save to file =====

if ~isempty(printpath)
    datnodat = datfile(1:end-4);
    print('-r300',[printpath 'calk_ptl_' datnodat],'-dpng')
end %if
    

end %function calk_ptl

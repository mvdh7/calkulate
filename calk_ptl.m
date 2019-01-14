%% Import data

% Initialise Python
% calk_initpy('//anaconda/envs/spritzer/bin/python') % Mac
calk_initpy( ... Windows
    'C:\Users\yau17reu\anaconda\Anaconda3\envs\spritzer\pythonw.exe')
% calk = py.importlib.import_module('calkulate');

% Settings
datfile = 'datfiles/0-0  0  (0)CRM-144-0435-4.dat';
Vsamp = 100;
Cacid = 0.1;
[CT,AT_cert,S,PT,SiT] = dicksonCRM(144);
CT      = CT      * 1e-6;
AT_cert = AT_cert * 1e-6;
PT      = PT      * 1e-6;
SiT     = SiT     * 1e-6;
burette_cx = 1;
Tk_force = [];

% Get initial Gran guesses
clear t
[t.Macid,t.EMF,t.Tk,Msamp, t.F1g,t.Lg,t.EMF0g, ATg,EMF0g,t.pHg,t.L] ...
    = calk_guessGran(datfile,Vsamp,Cacid,S);
t = struct2table(t);

t.mu = Msamp ./ (Msamp + t.Macid);

% Solve for TA fully
[AT,EMF0,AT_RMS,AT_Npts] = calk_VINDTA(datfile,Vsamp,Cacid, ...
    S,CT,PT,SiT,burette_cx,Tk_force);

[t.H,t.pH] = calk_EMF2H(t.EMF,EMF0,t.Tk);

[t.simAT,components] = calk_simAT(t.Macid,t.Tk,t.H,Msamp,S,CT,PT,SiT);
cvars = fieldnames(components);
t = [t struct2table(components)];

t.estAT = (t.simAT + t.Macid*Cacid ./ (t.Macid + Msamp)) ./ t.mu;

% ----------------------------------------------------- Plot the lot! -----

fvars = {'EMF' 'estAT' 'F1g' [] 'EMF0g'};
flabels = {'EMF / mV' [] 'F1' 'estAT' 'EMF0 guess / mV'};

fxlim = minmax(1e3*t.Macid');

figure(1); clf

for V = 1:numel(fvars)
if ~isempty(fvars{V})
    
subplot(3,2,V); hold on

    xlim(fxlim)

    scatter(1e3*t.Macid(t.Lg),t.(fvars{V})(t.Lg),'c','filled')
    scatter(1e3*t.Macid,t.(fvars{V}),'k')
    
    if ismember(fvars{V},'EMF')
        scatter(1e3*t.Macid(t.L),t.(fvars{V})(t.L),'r')
    end %if
    
    switch fvars{V}
        
        case 'EMF'
            
            ylim(minmax(t.EMF')+30*[-1 1])
            set(gca, 'ytick',0:100:800)
        
        case 'F1g'
            % plot regression line
            % plot ATg
        
        case 'EMF0g' 
            plot(fxlim,EMF0g*[1 1],'b')
            
        case 'estAT'
            plot(fxlim,AT*[1 1])
            
    end %switch

    xlabel('Acid mass / g')
    ylabel(flabels{V})
    
    set(gca, 'box','on', 'tickdir','out', 'xcolor','k', 'ycolor','k', ...
        'fontname','arial', 'fontsize',8, 'xtick',0:10)

end %if
end %for V

subplot(3,2,[4 6]); hold on
        
    plot(1e3*t.Macid,t.estAT)

    for V = 1:numel(cvars)

        plot(1e3*t.Macid,abs(t.(cvars{V})))

    end %for V

    legend(['estAT'; cvars], 'location','eastoutside')

    set(gca, 'yscale','log')

    xlim(fxlim)

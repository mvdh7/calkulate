%% Import data

% Initialise Python
calk_initpy('//anaconda/envs/spritzer/bin/python')
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

% Solve for TA fully
[AT,EMF0,AT_RMS,AT_Npts] = calk_VINDTA(datfile,Vsamp,Cacid, ...
    S,CT,PT,SiT,burette_cx,Tk_force);

[t.H,t.pH] = calk_EMF2H(t.EMF,EMF0,t.Tk);

% Calculate simulated TA
test = cell(py.calkulate.VINDTA.simAT(py.numpy.array(t.Macid'), ...
    py.numpy.array(t.Tk'),py.numpy.array(t.H'),Msamp,S,CT,PT,SiT));
t.AT = double(py.array.array('d',test{1}))';

% ----------------------------------------------------- Plot the lot! -----

fvars = {'EMF' 'F1g' 'EMF0g' 'AT'};
flabels = {'EMF / mV' 'F1' 'EMF0 guess / mV' 'AT'};

fxlim = minmax(t.Macid');

figure(1); clf

for V = 1:numel(fvars)

subplot(2,2,V); hold on

    scatter(t.Macid(t.Lg),t.(fvars{V})(t.Lg),'c','filled')
    scatter(t.Macid,t.(fvars{V}),'k')
    
    if ismember(fvars{V},'EMF')
        scatter(t.Macid(t.L),t.(fvars{V})(t.L),'r')
    end %if
    
    switch fvars{V}
        
        case 'F1g'
            % plot regression line
            % plot ATg
        
        case 'EMF0g' 
            plot(fxlim,EMF0g*[1 1],'b')
            
    end %switch

    xlim(fxlim)

    xlabel('Acid mass / g')
    ylabel(flabels{V})
    
end %for V
    
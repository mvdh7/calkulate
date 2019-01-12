%% Put your data here ... (change stuff on right hand side of the "=")

datfilepath = 'path/to/datfiles/'; % Path to .dat files

% CRMs
c_datfile = crmtable.datfile; % .dat filename from VINDTA
c_batch   = crmtable.batch;   % Dickson CRM batch number
c_bottle  = crmtable.bottle;  % Dickson CRM bottle number
c_acid    = crmtable.acid;    % Acid ID

% Samples
s_datfile = sampletable.datfile; % .dat filename from VINDTA
s_DIC     = sampletable.dic;     % Dissolved inorganic carbon / mol/kg
s_Phos    = sampletable.phos;    % Phosphate                  / mol/kg
s_Si      = sampletable.si;      % Silicate                   / mol/kg
s_Sal     = sampletable.s;       % Practical salinity
s_acid    = sampletable.acid;    % Acid ID

%% ... then run this (don't modify below here!)

% Test CRMs

for F = 1:numel(c_datfile)

    fid = fopen([datfilepath c_datfile{F}]);

    if fid > -1
        fclose(fid);
    else
        disp(['CRM .dat file ' num2str(F) ' not found: ' c_datfile{F}])
    end %if

end %for F
    

% Test samples

for F = 1:numel(s_datfile)

    fid = fopen(s_datfile{F});

    if fid > -1
        fclose(fid);
    else
        disp(['Sample .dat file ' num2str(F) ' not found: ' c_datfile{F}])
    end %if

end %for F

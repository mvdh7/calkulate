% Get list of datfiles
datlist = dir('datfiles/*.dat');

clear crm
crm.datfiles = cell(size(datlist));
crm = struct2table(crm);

for F = 1:height(crm)
    crm.datfiles{F} = datlist(F).name;
end %for F
clear datlist

% Identify batch numbers etc.
re_datfiles = regexp(crm.datfiles, ...
    'CRM-(\d{3})-(\d{4})-(\d)\.dat','tokens');

crm = addnancols(crm,{'batch' 'bottle' 'rep'});
for F = 1:height(crm)
    
    crm.batch (F) = str2double(re_datfiles{F}{1}{1});
    crm.bottle(F) = str2double(re_datfiles{F}{1}{2});
    crm.rep   (F) = str2double(re_datfiles{F}{1}{3});
    
end %for F

% Get certified values
[crm.dic_cert,crm.ta_cert,crm.s,crm.phos,crm.si] = dicksonCRM(crm.batch);

% Initialise Python
initpy('C:\Users\yau17reu\anaconda\Anaconda3\envs\calkenv\pythonw.exe')

% Get acid concentration
Vsamp = 100; % ml

crm = addnancols(crm,{'Cacid' 'ta' 'e0' 'ta_rms' 'npts'});
tic
for F = 1:height(crm)
    
    crm.Cacid(F) = calk_VINDTA_CRM(['datfiles/' crm.datfiles{F}], ...
        Vsamp,crm.ta_cert(F)*1e-6,crm.s(F),crm.dic_cert(F)*1e-6, ...
        crm.phos(F)*1e-6,crm.si(F)*1e-6,1,[]);
    
end %for F
toc

Cacid = median(crm.Cacid);

% Calibrate samples
tic
for F = 1:height(crm)
    
%     Fsamp = py.calkulate.VINDTA.MPH(['datfiles/' crm.datfiles{F}], ...
%         Vsamp,Cacid,crm.s(F),crm.dic_cert(F)*1e-6, ...
%         crm.phos(F)*1e-6,crm.si(F)*1e-6,1);
%     
%     Fsamp = numpy2double(Fsamp{'x'});
    
    [crm.ta(F),crm.e0(F),crm.ta_rms(F),crm.npts(F)] ...
        = calk_VINDTA(['datfiles/' crm.datfiles{F}], ...
        Vsamp,Cacid,crm.s(F),crm.dic_cert(F)*1e-6, ...
        crm.phos(F)*1e-6,crm.si(F)*1e-6,1,[]);
    
end %for F
toc

% Convert units
crm.ta     = crm.ta     * 1e6;
crm.ta_rms = crm.ta_rms * 1e6;

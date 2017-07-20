function [acidmolar,g_acidmolar,g_ta] = CalkCRM(datfile, ...
    acidmolar_approx,acidrho,pipvol,CRMbatch,burette_cxn,sc,tforce)
%CalkCRM Gets best HCl molarity for a total alkalinity certified reference
%  material (CRM) titration.
% 
% This works for one titration at a time. How to use the output is left to
%  the user's discretion (e.g. average by acid batch, or analysis session).
% 
% Instead of using DicksonCRM function, you can set <cert_dic>, <cert_ta>,
%  <cert_s> and <cert_phos> in script to certified values. The function
%  returns NaN if CRM batch data is not in DicksonCRM - just add it in (at
%  the end of this script)!
% 
% Output <acidmolar> becomes input <acidmolar> in function calkulate.
% 
% Inputs:  datfile = VINDTA .dat file name for CRM
% acidmolar_approx = Approximate acid molarity / mol/l
%          acidrho = Acid density / kg/l
%           pipvol = TA pipette volume / ml
%         CRMbatch = Dickson CRM batch number
% NB: requires function dicksonCRM, and that function must be contain 
%     information for the batch in question, else a NaN will be returned.
%      burette_cxn = Titrino pipette volume correction factor
%               sc = structure specifying all constants' codes (can be left
%                    empty)
%           tforce = Leave empty to use titration file temperature, else
%                    provide temp. to over-ride with / degC
% 
% NB: the same value for "acidrho", "pipvol", "sc" and "burette_cxn" must 
%  then be used by Calkulate. If any of them are unknown, use estimated
%  values; the calculated acid molarity should adjust itself to take these 
%  things into account too. In that case, it is no longer really an acid
%  molarity any more, but a generic correction factor taking all of these
%  things into account. If they are unknown, it is essential that these
%   volumes etc. do not change between the CRM titration and the subsequent
%  sample titrations.
% 
% Written by Matthew P. Humphreys, 2015-01-21. Last updated 2016-07-07.
% === Please cite as: ===
% Humphreys, M.P. (2015), "Calculating seawater total alkalinity from
%  open-cell titration data using a modified Gran plot technique", in
%  "Measurements and Concepts in Marine Carbonate Chemistry", PhD thesis,
%  Ocean and Earth Science, University of Southampton, UK, pp. 25-44.

%% Test inputs
% datfile = vcrm_datfile{C};
% acidmolar_approx = 0.1;
% acidrho = vindta_v38_acidrho;
% pipvol = vindta_v38_pipvol;
% CRMbatch = vcrm_batch(C);
% burette_cxn = vindta_v38_tit5vol/5;
% sc = [];
% tforce = [];

%% Certified TA
[cert_dic,cert_ta,cert_s,cert_phos] = DicksonCRM(CRMbatch);

if ~isnan(cert_ta)

final_value = false;
g_acidmolar = linspace(0.8*acidmolar_approx,1.2*acidmolar_approx,5);

acidmolar_prev = NaN;

while ~final_value
    g_ta = NaN(size(g_acidmolar));
    
for C = 1:length(g_acidmolar)
    g_ta(C) = Calkulate(datfile,g_acidmolar(C),acidrho, ...
        pipvol,cert_s,cert_dic*1e-6,cert_phos*1e-6,burette_cxn,sc, ...
        tforce,true);
end %for C

amfit1 = regstats(g_acidmolar,g_ta,'quadratic',{'beta' 'yhat'});
first_am = polyval(flipud(amfit1.beta)',cert_ta);

acidmolar = abs(first_am);

% disp([acidmolar_prev acidmolar_approx acidmolar]);

if acidmolar < max(g_acidmolar) && acidmolar > min(g_acidmolar) ...
	&& (abs(diff([acidmolar acidmolar_approx])) ...
        < acidmolar_approx * 1e-5 ...
        || abs(diff([acidmolar acidmolar_prev])) ...
        < acidmolar_prev * 1e-5)
    final_value = true;
    acidmolar = mean([acidmolar acidmolar_approx]);
else
    acidmolar_prev = acidmolar_approx;
    acidmolar_approx = acidmolar;
    g_acidmolar = linspace(0.99*acidmolar_approx, ...
        1.01*acidmolar_approx,5);
end %if else
    
end %while

else
    acidmolar = NaN;
end %if else
    
end %function CalkCRM

function [dic,ta,s,phos,si,no2,no3,batchCert] = DicksonCRM(batches)
%dicksonCRM Gets certified values for carbonate chemistry CRMs
%  from AG Dickson (Scripps Institution of Oceanography, USA)
% Input <batches> must be numerical and a single value or column vector
% Outputs: <dic> Dissolved Inorganic Carbon / micromol/kg
%           <ta> Total Alkalinity           / micromol/kg
%            <s> Salinity
%         <phos> Phosphate                  / micromol/kg
%           <si> Silicate                   / micromol/kg
%          <no2> Nitrite                    / micromol/kg
%          <no3> Nitrate                    / micromol/kg
%    <batchCert> All of above data in one matrix
% Source: http://cdiac.esd.ornl.gov/oceans/Dickson_CRM/batches.html
% Written by Matthew P. Humphreys, last updated 2015-06-05.

%% Certified batch values - incomplete list!
%       Batch    DIC       TA     Salinity  Phos. Silic.  NO2   NO3
cert = [ 104 , 2020.10 , 2222.61 , 33.300 , 0.77 , 4.4 , 0.00 , 6.70 ;
         105 , 2036.72 , 2235.16 , 33.504 , 0.62 , 4.9 , 0.01 , 3.00 ;
         107 , 1995.12 , 2205.17 , 32.985 , 0.37 , 2.7 , 0.00 , 0.57 ;
         109 , 2026.33 , 2224.26 , 33.328 , 0.40 , 2.8 , 0.00 , 2.52 ;
         114 , 2000.93 , 2217.91 , 33.208 , 0.36 , 2.4 , 0.00 , 0.97 ;
         117 , 2009.99 , 2239.18 , 33.503 , 0.37 , 1.7 , 0.00 , 1.14 ;
         119 , 2014.78 , 2221.24 , 33.263 , 0.34 , 2.8 , 0.00 , 1.21 ;
         120 , 2002.61 , 2208.43 , 33.072 , 0.31 , 2.8 , 0.00 , 0.23 ;
         121 , 2039.26 , 2225.01 , 33.346 , 0.77 , 8.7 , 0.00 , 6.10 ;
         123 , 2022.04 , 2225.21 , 33.384 , 0.48 , 4.8 , 0.00 , 2.70 ;
         128 , 2013.54 , 2240.28 , 33.442 , 0.30 , 3.7 , 0.00 , 0.22 ;
         133 , 2021.12 , 2224.37 , 33.341 , 0.40 , 2.9 , 0.00 , 0.70 ;
         134 , 2026.91 , 2236.51 , 33.651 , 0.42 , 3.3 , 0.00 , 2.90 ;
         135 , 2036.96 , 2226.33 , 33.561 , 0.55 , 2.3 , 0.00 , 2.50 ;
         136 , 2021.15 , 2246.74 , 33.678 , 0.33 , 1.4 , 0.00 , 0.11 ;
         138 , 2033.98 , 2221.75 , 33.398 , 0.58 , 1.5 , 0.00 , 1.69 ;
         141 , 2033.26 , 2234.07 , 33.417 , 0.48 , 2.7 , 0.02 , 2.39 ;
         142 , 2038.07 , 2227.59 , 33.389 , 0.29 , 3.3 , 0.02 , 4.63 ;
         144 , 2031.53 , 2238.60 , 33.571 , 0.31 , 2.5 , 0.00 , 0.85 ;
         146 , 2002.92 , 2214.11 , 33.122 , 0.34 , 1.8 , 0.01 , 0.22 ;
         148 , 2024.12 , 2241.13 , 33.624 , 0.26 , 2.4 , 0.00 , 0.35 ;
         151 , 2033.83 , 2225.56 , 33.345 , 0.56 , 3.5 , 0.00 , 1.53 ;
         152 , 2020.88 , 2216.94 , 33.371 , 0.51 , 5.7 , 0.01 , 1.30 ;
         154 , 2037.68 , 2224.30 , 33.347 , 0.61 , 4.4 , 0.00 , 2.30];

%% Get certified values for input batch(es)
batchCert = NaN(size(batches,1),7);
for U = 1:length(batches)
    if any(cert(batches(U) == cert(:,1)))
        batchCert(U,:) = cert(batches(U) == cert(:,1),2:end);
    else
        if ~isnan(batches(U))
            disp(['WARNING! Data for Dickson CRM batch ' ...
                num2str(batches(U)) ...
                ' is not yet stored in this function']);
        end %if
    end %if elseif
end %for U

%% Split up result for output
dic  = batchCert(:,1);
ta   = batchCert(:,2);
s    = batchCert(:,3);
phos = batchCert(:,4);
si   = batchCert(:,5);
no2  = batchCert(:,6);
no3  = batchCert(:,7);

end %function DicksonCRM
function [AT,components] = calk_simAT(Macid,Tk,H,Msamp,S,CT,PT,SiT)

% Calculate simulated TA
simAT = cell(py.calkulate.VINDTA.simAT(py.numpy.array(Macid(:)'), ...
    py.numpy.array(Tk(:)'),py.numpy.array(H(:)'),Msamp,S,CT,PT,SiT));
AT = double(py.array.array('d',simAT{1}))';

% Get TA components
simComp = cell(simAT{2});

cvars = {'bicarb' 'carb' 'B4' 'OH' [] 'HSO4' 'HF' 'P0' 'P2' 'P3' 'SiOOH3'};
for V = 1:numel(cvars)
    if ~isempty(cvars{V})
        components.(cvars{V}) = double(py.array.array('d',simComp{V}))';
    end %if
end %for V

end %function calk_simAT

function [H,pH] = calk_simH(Macid,Tk,Msamp,Cacid,S,AT,CT,PT,SiT)

H = py.calkulate.VINDTA.simH(py.numpy.array(Macid(:)'), ...
    py.numpy.array(Tk(:)'),Msamp,Cacid,S,AT,CT,PT,SiT);

H = double(py.array.array('d',H))';

pH = -log10(H);

end %function calk_simH

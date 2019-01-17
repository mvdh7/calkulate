function [H,pH] = calk_EMF2H(EMF,EMF0,Tk)

H = py.calkulate.solve.EMF2H(py.numpy.array(EMF(:)'),EMF0, ...
    py.numpy.array(Tk(:)'));

H = double(py.array.array('d',H))';

pH = -log10(H);

end %function calk_EMF2H

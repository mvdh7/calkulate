function [Vacid,EMF,Tk] = calk_getdat(datfile)

data = cell(py.calkulate.gettit.VINDTA(datfile));

Vacid = double(py.array.array('d',data{1}))';
EMF   = double(py.array.array('d',data{2}))';
Tk    = double(py.array.array('d',data{3}))';

end %function calk_getdat

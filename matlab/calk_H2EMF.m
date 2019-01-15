function EMF = calk_H2EMF(H,EMF0,Tk)

EMF = py.calkulate.solve.H2EMF(py.numpy.array(H(:)'),EMF0, ...
    py.numpy.array(Tk(:)'));

EMF = double(py.array.array('d',EMF))';

end %function calk_H2EMF

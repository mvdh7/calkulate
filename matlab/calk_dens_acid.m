function rho_acid = calk_dens_acid(Tk)

rho_acid = py.calkulate.dens.acid(py.numpy.array(Tk(:)'));

rho_acid = double(py.array.array('d',rho_acid))';

end %function calk_dens_acid

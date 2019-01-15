function rho_sw = calk_dens_sw(Tk,S)

rho_sw = py.calkulate.dens.sw(py.numpy.array(Tk(:)'), ...
    py.numpy.array(S(:)'));

rho_sw = double(py.array.array('d',rho_sw))';

end %function calk_dens_sw

function [Macid,EMF,Tk,Msamp, F1g,Lg,EMF0gvec, ATg,EMF0g,pHg,L] ...
    = calk_guessGran(datfile,Vsamp,Cacid,S)

guessGran = cell(py.calkulate.VINDTA.guessGran(datfile,Vsamp,Cacid,S));

Macid    =  double(py.array.array('d',guessGran{ 1}))' * 1e3;
EMF      =  double(py.array.array('d',guessGran{ 2}))';
Tk       =  double(py.array.array('d',guessGran{ 3}))';
Msamp    =                            guessGran{ 4}   ;
F1g      =  double(py.array.array('d',guessGran{ 5}))';
Lg       = logical(py.array.array('d',guessGran{ 6}))';
EMF0gvec =  double(py.array.array('d',guessGran{ 7}))';
ATg      =                            guessGran{ 8}   ;
EMF0g    =                            guessGran{ 9}   ;
pHg      =  double(py.array.array('d',guessGran{10}))';
L        = logical(py.array.array('d',guessGran{11}))';

end %function calk_guessGran

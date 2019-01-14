function calk_nolegend(handle)
%calk_nolegend Removes legend entry for plot specified by <handle>
set(get(get(handle,'annotation'),'legendinformation'), ...
        'icondisplaystyle','off');
end %function calk_nolegend

function m = temporal_moment(t,P_y,i)
    m = trapz(t,(t.^i.*P_y))
end
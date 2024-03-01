using Plots

function fscss_signal(sf::Int64,bw::Int64,m::Int64;ovs=1)
    M=2^sf #number of chips
    tc=1/bw # chip length [sec]
    ts=tc*M # symbol length [sec]
    M=M*ovs # if you want to see signal form that increase instantaneous frequency, you should set "ovs=2".  
    x=zeros(ComplexF64,M) # fscc (LoRa) signal on equivalent low-pass system

    for k in 1:M
        x[k]=exp(im*2pi*mod((m*ovs+(k-1))*tc/ovs,M*tc)^2*bw/(2*M*tc))
    end

    axis_x=range(0,ts,length=M)
    p1=plot(axis_x,real(x),
        legend=false,
        xgrid=false,
        ygrid=false,
        fmt=:svg)
    plot!(axis_x,imag(x),
        legend=false,
        xgrid=false,
        ygrid=false,)
    savefig(p1,"fscss_signal_m$(m).svg")
    display(p1)
    println("Press Enter to exit this program (or function)")
    readline(stdin)
end

fscss_signal(7,125*10^3,0;ovs=1)
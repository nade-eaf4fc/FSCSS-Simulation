using Plots

function fscss_signal(sf::Int64,bw::Int64,m::Int64;ovs=1)
    M=2^sf
    tc=1/bw
    ts=tc*M
    M=M*ovs
    x=zeros(ComplexF64,M)
    x0=zeros(ComplexF64,M)

    for k in 1:M
        x[k]=exp(im*2pi*mod((m+ovs*(k-1))*tc/ovs,M*tc)^2*bw/(2*M*tc))
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
    savefig(p1,"fscss_signal.svg")
    display(p1)

    
    readline(stdin)
end

fscss_signal(7,125*10^3,0)
using Plots
using FFTW
gr()
function dft(y::Array,M::Int64;ovs=1)
    Y=zeros(ComplexF64,M)#周波数領域
    for n in 0:M-1
        for k in 0:M-1 #配列は1から始まることに注意
            Y[n+1] += y[k+1]*exp(-im*2pi*n*k/M)
        end
    end
    return Y
end

function fscss_demodu_dft(sf::Int64,bw::Int64;ovs=1)
    M=2^sf #number of chips
    tc=1/bw # chip length [sec]
    ts=tc*M # symbol length [sec]
    M=M*ovs # if you want to see signal form that increase instantaneous frequency, you should set "ovs=2".  
    x=zeros(ComplexF64,M) # fscc (LoRa) signal on equivalent low-pass system
    Y=zeros(ComplexF64,M)
    d=zeros(ComplexF64,M)
    anime=Animation()


    for m in 0:2^sf-1

        for k in 1:M
            x[k]=exp(im*2pi*mod((m*ovs+(k-1))*tc/ovs,M*tc)^2*bw/(2*ts))
        end

        for k in 1:M
            d[k]=x[k]*exp(-im*2pi*((k-1)*tc/ovs)^2*bw/(2*ts))
        end

        Y=dft(d,M,ovs=ovs)
        Y_abs=1/(M)*abs.(Y)
        axis_x=range(0,M-1,length=M)
        p1=plot(axis_x,Y_abs,
        xgrid=false,
        ygrid=false,
        lc="#000000",
        xlims=(0,M),
        ylims=(0,1),
        legend=false)
        
        frame(anime,p1)

    end

    gif(anime,"dftout_sf$(sf)_ovs$(ovs).gif";fps=20)
    println("Press Enter to exit this program (or function)")
    readline(stdin)
end

fscss_demodu_dft(7,125*10^3;ovs=4)
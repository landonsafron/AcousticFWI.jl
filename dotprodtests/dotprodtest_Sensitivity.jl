using AcousticFWI,Seismic,PyPlot

function main()

    ns = 4
    ng = 250
    nz = 150
    nx = 250
    nt = 1024
    nf = 1024
    dz = 0.01
    dx = 0.01
    dt = 0.002
    f0 = 10.
    fmin = 0.5 
    fmax = 30.
    ext = 50
    atten_max = 2.5

    vp = 2.1*ones(nz,nx)
    vp[61:100,:] = 3.0
    vp[101:end,:] = 2.7
    
    vp0 = 2.5*ones(vp)

    ot = zeros(ns)
    isz = 3*ones(Int,ns)
    isx = round(Int,collect(linspace(2,nx-1,ns)))
    igz = 3*ones(Int,ng)
    igx = round(Int,collect(linspace(1,nx,ng)))

    wav = Ricker(f0=f0,dt=dt)
    d = zeros(nt,ng,ns)
    for ishot = 1:ns
        u = Helmholtz(isz[ishot],isx[ishot],ot[ishot],vp,wav,fmin,fmax,nf,nt,dz,dx,dt,ext,atten_max)
        d[:,:,ishot] = u[:,3,igx]
    end


    R = Restriction(nz,nx,ext,igz,igx)

end
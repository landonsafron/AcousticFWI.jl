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

    vp0 = 2.5*ones(vp)

    ## Perform 1 iteration of FWI to check if the Sensitivity operator passes the dot-product test
    
    nz = size(vp0,1) + 2*ext
    nx = size(vp0,2) + 2*ext

    L = Laplacian(nz,nx,dz,dx)
    M = MassMatrix(vp0,nz,nx,ext)
    A = Attenuation(nz,nx,ext,atten_max)
    R = Restriction(nz,nx,ext,igz,igx)

    D = fft([d ; zeros(nf-nt,ng,ns)],1)
    WAV = fft([wav;zeros(nf-length(wav))])
    fs = 1/dt
    df = fs/nf
    faxis = fftshift(-fs/2:df:fs/2-df)
    waxis = 2*pi*faxis

    _,iwmin = findmin(abs(fmin-faxis))
    _,iwmax = findmin(abs(fmax-faxis))

    param = Dict(:R=>R,:A=>A)

    iw = 15
    w = waxis[iw]
    s = zeros(eltype(D),nz*nx,ns)
    for ishot = 1:ns
        s[:,ishot] = Source(isz[ishot],isx[ishot],ot[ishot],WAV,waxis,w,nz,nx,ext)
    end
    H = lufact(L + w^2*M*A)
    U = zeros(eltype(D),nz*nx,ns)
    r = zeros(eltype(D),ng*ns)
    for ishot = 1:ns
        U[:,ishot] = H\s[:,ishot]
        r[(ishot-1)*ng+1:ishot*ng] = R*U[:,ishot] - D[iw,:,ishot]
    end
    param[:H] = H
    param[:U] = U
    param[:w] = w

    m = rand(nz*nx) + im*rand(nz*nx)
    r = rand(ns*ng) + im*rand(ns*ng)

    R = SensitivityMultiShot(m,false;param...)
    M = SensitivityMultiShot(r,true;param...)

    dotprodtest = dot(m[:],M[:]) - dot(r[:],R[:])
    println("This number should be close to zero: ",dotprodtest)

end

main()
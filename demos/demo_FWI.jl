using AcousticFWI,Seismic,PyPlot

#function main()

    ns = 5  # ns = 10
    ng = 250
    nz = 150
    nx = 250
    nt = 1024
    nf = 1024
    dz = 0.01
    dx = 0.01
    dt = 0.002
    f0 = 10.
    fmin = 1. 
    fmax = 30.
    ext = 50
    atten_max = 2.5
    GNiter = 4
    CGiter = 20

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
    vp_est = FWI(vp0,d,wav,isz,isx,igz,igx,ot,fmin,fmax,nf,dz,dx,dt,ext,atten_max,GCiter,CGiter)

    subplot(3,3,1) ; SeisPlot(u[10,:,:],pclip=100,cmap="gray",fignum=1)
    subplot(3,3,2) ; SeisPlot(u[70,:,:],pclip=100,cmap="gray",fignum=1)
    subplot(3,3,3) ; SeisPlot(u[130,:,:],pclip=100,cmap="gray",fignum=1)
    subplot(3,3,4) ; SeisPlot(u[190,:,:],pclip=100,cmap="gray",fignum=1)
    subplot(3,3,5) ; SeisPlot(u[250,:,:],pclip=100,cmap="gray",fignum=1)
    subplot(3,3,6) ; SeisPlot(u[310,:,:],pclip=100,cmap="gray",fignum=1)
    subplot(3,3,7) ; SeisPlot(u[370,:,:],pclip=100,cmap="gray",fignum=1)
    subplot(3,3,8) ; SeisPlot(u[430,:,:],pclip=100,cmap="gray",fignum=1)
    subplot(3,3,9) ; SeisPlot(u[490,:,:],pclip=100,cmap="gray",fignum=1)

    SeisPlot(d,pclip=98,cmap="gray",fignum=2)

    figure(3) ; vmin = min(minimum(vp),minimum(vp0),minimum(vp_est))
              ; vmax = max(maximum(vp),maximum(vp0),maximum(vp_est))
              ; subplot(311) ; imshow(vp,vmin=vmin,vmax=vmax,cmap="YlOrBr",extent=[0.,(nx-1)*dx,(nz-1)*dz,0.])
              ; subplot(312) ; imshow(vp0,vmin=vmin,vmax=vmax,cmap="YlOrBr",extent=[0.,(nx-1)*dx,(nz-1)*dz,0.])
              ; subplot(313) ; imshow(vp_est,vmin=vmin,vmax=vmax,cmap="YlOrBr",extent=[0.,(nx-1)*dx,(nz-1)*dz,0.])

#end
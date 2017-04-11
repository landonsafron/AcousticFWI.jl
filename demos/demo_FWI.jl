using AcousticFWI,Seismic,SeismicImaging,PyPlot

    nz = 150
    nx = 250
    nt = 500
    nf = 1024
    dz = 2.
    dx = 2.
    dt = 0.002
    f0 = 45.
    fmin = 2. 
    fmax = 100.
    ext = 50
    atten_max = 2.
    alpha = 0.01
    maxiter = 2

    vp = 2000.*ones(nz,nx)
    vp[101:end,:] = 3000.

    isz = [3,3]
    isx = [81,167]
    igz = 3*ones(Int,nx)
    igx = [1:nx;]
    ot = [0.,30*dt]

    wav = Ricker(f0=f0,dt=dt)
    u = HelmholtzSolver(isz,isx,ot,vp,wav,fmin,fmax,nf,nt,dz,dx,dt,ext,atten_max)
    d = u[:,3,igx]

    vp0 = smooth2d(vp,50,50)
    vp_est = FWI(vp0,d,wav,isz,isx,igz,igx,ot,fmin,fmax,nf,dz,dx,dt,ext,atten_max,maxiter)

    subplot(4,4,1) ; SeisPlot(u[20,:,:],pclip=100,cmap="gray",fignum=1)
    subplot(4,4,2) ; SeisPlot(u[40,:,:],pclip=100,cmap="gray",fignum=1)
    subplot(4,4,3) ; SeisPlot(u[60,:,:],pclip=100,cmap="gray",fignum=1)
    subplot(4,4,4) ; SeisPlot(u[80,:,:],pclip=100,cmap="gray",fignum=1)
    subplot(4,4,5) ; SeisPlot(u[100,:,:],pclip=100,cmap="gray",fignum=1)
    subplot(4,4,6) ; SeisPlot(u[120,:,:],pclip=100,cmap="gray",fignum=1)
    subplot(4,4,7) ; SeisPlot(u[140,:,:],pclip=100,cmap="gray",fignum=1)
    subplot(4,4,8) ; SeisPlot(u[160,:,:],pclip=100,cmap="gray",fignum=1)
    subplot(4,4,9) ; SeisPlot(u[180,:,:],pclip=100,cmap="gray",fignum=1)
    subplot(4,4,10) ; SeisPlot(u[200,:,:],pclip=100,cmap="gray",fignum=1)
    subplot(4,4,11) ; SeisPlot(u[220,:,:],pclip=100,cmap="gray",fignum=1)
    subplot(4,4,12) ; SeisPlot(u[240,:,:],pclip=100,cmap="gray",fignum=1)
    subplot(4,4,13) ; SeisPlot(u[260,:,:],pclip=100,cmap="gray",fignum=1)
    subplot(4,4,14) ; SeisPlot(u[280,:,:],pclip=100,cmap="gray",fignum=1)
    subplot(4,4,15) ; SeisPlot(u[300,:,:],pclip=100,cmap="gray",fignum=1)
    subplot(4,4,16) ; SeisPlot(u[320,:,:],pclip=100,cmap="gray",fignum=1)

    SeisPlot(d,pclip=98,cmap="gray",fignum=2)

    figure(3) ; vmin = min(minimum(vp),minimum(vp0),minimum(vp_est))
              ; vmax = max(maximum(vp),maximum(vp0),maximum(vp_est))
              ; subplot(311) ; imshow(vp,vmin=vmin,vmax=vmax,cmap="YlOrBr",extent=[0.,(nx-1)*dx,(nz-1)*dz,0.])
              ; subplot(312) ; imshow(vp0,vmin=vmin,vmax=vmax,cmap="YlOrBr",extent=[0.,(nx-1)*dx,(nz-1)*dz,0.])
              ; subplot(313) ; imshow(vp_est,vmin=vmin,vmax=vmax,cmap="YlOrBr",extent=[0.,(nx-1)*dx,(nz-1)*dz,0.])

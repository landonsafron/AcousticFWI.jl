using AcousticFWI,Seismic,SeismicImaging,PyPlot

#function main()

    nz = 150
    nx = 250
    nt = 512
    nf = 512
    dz = 2.
    dx = 2.
    dt = 0.002
    f0 = 45.
    fmin = 2. 
    fmax = 120.
    ext = 50
    atten_max = 2.
    alpha = 1.0e-8
    maxiter = 2

    vp = 2000.*ones(nz,nx)
    vp[101:end,:] = 3000.

    ot = [0.,0.,0.,0.]
    isz = [3,3,3,3]
    isx = [50,100,150,200]
    igz = 3*ones(Int,nx)
    igx = [1:nx;]

    wav = Ricker(f0=f0,dt=dt)
    u = HelmholtzSolver(isz,isx,ot,vp,wav,fmin,fmax,nf,nt,dz,dx,dt,ext,atten_max)
    d = u[:,3,igx]

    vp0 = smooth2d(vp,20,20)
    vp_est = FWI(vp0,d,wav,isz,isx,igz,igx,ot,fmin,fmax,nf,dz,dx,dt,ext,atten_max,maxiter)

    subplot(3,3,1) ; SeisPlot(u[20,:,:],pclip=100,cmap="gray",fignum=1)
    subplot(3,3,2) ; SeisPlot(u[40,:,:],pclip=100,cmap="gray",fignum=1)
    subplot(3,3,3) ; SeisPlot(u[60,:,:],pclip=100,cmap="gray",fignum=1)
    subplot(3,3,4) ; SeisPlot(u[80,:,:],pclip=100,cmap="gray",fignum=1)
    subplot(3,3,5) ; SeisPlot(u[100,:,:],pclip=100,cmap="gray",fignum=1)
    subplot(3,3,6) ; SeisPlot(u[120,:,:],pclip=100,cmap="gray",fignum=1)
    subplot(3,3,7) ; SeisPlot(u[140,:,:],pclip=100,cmap="gray",fignum=1)
    subplot(3,3,8) ; SeisPlot(u[160,:,:],pclip=100,cmap="gray",fignum=1)
    subplot(3,3,9) ; SeisPlot(u[180,:,:],pclip=100,cmap="gray",fignum=1)

    SeisPlot(d,pclip=98,cmap="gray",fignum=2)

    figure(3) ; vmin = min(minimum(vp),minimum(vp0),minimum(vp_est))
              ; vmax = max(maximum(vp),maximum(vp0),maximum(vp_est))
              ; subplot(311) ; imshow(vp,vmin=vmin,vmax=vmax,cmap="YlOrBr",extent=[0.,(nx-1)*dx,(nz-1)*dz,0.])
              ; subplot(312) ; imshow(vp0,vmin=vmin,vmax=vmax,cmap="YlOrBr",extent=[0.,(nx-1)*dx,(nz-1)*dz,0.])
              ; subplot(313) ; imshow(vp_est,vmin=vmin,vmax=vmax,cmap="YlOrBr",extent=[0.,(nx-1)*dx,(nz-1)*dz,0.])

#end
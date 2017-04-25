using AcousticFWI,Seismic,PyPlot

function main()

    ns = 10
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
    GNiter = 10
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
    vp_est = FWI(vp0,d,wav,isz,isx,igz,igx,ot,fmin,fmax,nf,dz,dx,dt,ext,atten_max,GNiter,CGiter)

    figure(1) ; vmin = 1000*min(minimum(vp),minimum(vp0),minimum(vp_est))
              ; vmax = 1000*max(maximum(vp),maximum(vp0),maximum(vp_est))
              ; subplot(311) ; imshow(1000*vp,vmin=vmin,vmax=vmax,cmap="YlOrBr",extent=[0.,(nx-1)*1000*dx,(nz-1)*1000*dz,0.])
                             ; xlabel("x (m)")
                             ; ylabel("z (m)")
                             ; colorbar(aspect=15,ticks=[2100,2300,2500,2700,2900,3100])
                             ; ax = gca() 
                             ; ax[:text](-600.,-100., "a)", fontsize=15, fontweight="bold")
              ; subplot(312) ; imshow(1000*vp0,vmin=vmin,vmax=vmax,cmap="YlOrBr",extent=[0.,(nx-1)*1000*dx,(nz-1)*1000*dz,0.])
                             ; xlabel("x (m)")
                             ; ylabel("z (m)")
                             ; colorbar(aspect=15,ticks=[2100,2300,2500,2700,2900,3100])
                             ; ax = gca() 
                             ; ax[:text](-600.,-100., "b)", fontsize=15, fontweight="bold")
              ; subplot(313) ; imshow(1000*vp_est,vmin=vmin,vmax=vmax,cmap="YlOrBr",extent=[0.,(nx-1)*1000*dx,(nz-1)*1000*dz,0.])
                             ; xlabel("x (m)")
                             ; ylabel("z (m)")
                             ; colorbar(aspect=15,ticks=[2100,2300,2500,2700,2900,3100])
                             ; ax = gca() 
                             ; ax[:text](-600.,-100., "c)", fontsize=15, fontweight="bold")

end
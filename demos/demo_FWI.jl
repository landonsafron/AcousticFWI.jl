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
    GNiter = 8
    CGiter = 20

    vp = 2.1*ones(nz,nx)
    vp[61:100,:] = 3.0
    vp[101:end,:] = 2.5
    for (i,iz) in enumerate(61:100)
        vp[iz,1:105+i] = 2.7
    end

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

    vp0 = 2.3*ones(vp)
    vp_est,cost,f_inv = FWI(vp0,d,wav,isz,isx,igz,igx,ot,fmin,fmax,nf,dz,dx,dt,ext,atten_max,GNiter,CGiter)

    SeisPlot(d[:,:,round(Int,ns/2)],pclip=98,dx=1,dy=dt,xlabel="Trace",ylabel="t (s)",cmap="gray")
    figure(2) ; vmin = 1000*min(minimum(vp),minimum(vp0),minimum(vp_est))
              ; vmax = 1000*max(maximum(vp),maximum(vp0),maximum(vp_est))
              ; subplot(221) ; imshow(1000*vp,vmin=vmin,vmax=vmax,cmap="YlOrBr",extent=[0.,(nx-1)*1000*dx,(nz-1)*1000*dz,0.])
                             ; xlabel("x (m)")
                             ; ylabel("z (m)")
                             ; colorbar(aspect=15,ticks=[2100,2300,2500,2700,2900,3100])
                             ; ax = gca() 
                             ; ax[:text](-600.,-100., "a)", fontsize=15, fontweight="bold")
                             ; ax[:set_xticks]([0,500,1000,1500,2000])
                             ; ax[:set_yticks]([0,300,600,900,1200,1500])
              ; subplot(222) ; imshow(1000*vp0,vmin=vmin,vmax=vmax,cmap="YlOrBr",extent=[0.,(nx-1)*1000*dx,(nz-1)*1000*dz,0.])
                             ; xlabel("x (m)")
                             ; ylabel("z (m)")
                             ; colorbar(aspect=15,ticks=[2100,2300,2500,2700,2900,3100])
                             ; ax = gca() 
                             ; ax[:text](-600.,-100., "b)", fontsize=15, fontweight="bold")
                             ; ax[:set_xticks]([0,500,1000,1500,2000])
                             ; ax[:set_yticks]([0,300,600,900,1200,1500])
              ; subplot(223) ; imshow(1000*vp_est,vmin=vmin,vmax=vmax,cmap="YlOrBr",extent=[0.,(nx-1)*1000*dx,(nz-1)*1000*dz,0.])
                             ; xlabel("x (m)")
                             ; ylabel("z (m)")
                             ; colorbar(aspect=15,ticks=[2100,2300,2500,2700,2900,3100])
                             ; ax = gca() 
                             ; ax[:text](-600.,-100., "c)", fontsize=15, fontweight="bold")
                             ; ax[:set_xticks]([0,500,1000,1500,2000])
                             ; ax[:set_yticks]([0,300,600,900,1200,1500])
              ; subplot(224) ; plot(1:length(cost[1]),cost[1],color="black",linewidth=2)
                             ; xlabel("Iteration")
                             ; ylabel("Normalized Misfit")
                             ; ax = gca() 
                             ; ax[:text](-0.4,1.05, "d)", fontsize=15, fontweight="bold")
                             ; ax[:set_xticks]([1,2,3,4,5,6,7,8])
                             ; ax[:set_yticks]([0.2,0.4,0.6,0.8,1.0])

end

main()
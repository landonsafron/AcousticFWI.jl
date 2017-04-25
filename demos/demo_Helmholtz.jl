using AcousticFWI,Seismic,PyPlot

function main()

    nz = 150
    nx = 250
    nt = 1024
    nf = 1024
    dz = 0.01
    dx = 0.01
    dt = 0.002
    f0 = 10.
    fmin = 1.0 
    fmax = 30.
    ext = 50
    atten_max = 2.5

    vp = 2.1*ones(nz,nx)
    vp[61:100,:] = 3.0
    vp[101:end,:] = 2.7

    #isz = [3]
    #isx = [125]
    #igz = 3*ones(Int,nx)
    #igx = [1:nx;]
    #ot = [0.,30*dt]
    isz = [3,3]
    isx = [80,170]
    igz = 3*ones(Int,nx)
    igx = [1:nx;]
    ot = [0.,64*dt]

    wav = Ricker(f0=f0,dt=dt)
    u = Helmholtz(isz,isx,ot,vp,wav,fmin,fmax,nf,nt,dz,dx,dt,ext,atten_max)
    d = u[:,3,igx]

    figure(1) ; imshow(1000*vp,cmap="YlOrBr",extent=[0.,(nx-1)*1000*dx,(nz-1)*1000*dz,0.])
              ; xlabel("x (m)")
              ; ylabel("z (m)")
              ; colorbar(ticks=[2100,2300,2500,2700,2900])

    figure(2) ; subplots_adjust(wspace=0.5,hspace=0.5)
                subplot(3,3,1) ; SeisPlot(u[10,:,:],pclip=100,dx=1000*dx,dy=1000*dz,xticks=[0,500,1000,1500,2000],yticks=[0,300,600,900,1200],xlabel="x (m)",ylabel="z (m)",title="t = 0.02s",titlesize=14,labelsize=12,cmap="gray",fignum=2)
                subplot(3,3,2) ; SeisPlot(u[75,:,:],pclip=100,dx=1000*dx,dy=1000*dz,xticks=[0,500,1000,1500,2000],yticks=[0,300,600,900,1200],xlabel="x (m)",ylabel="z (m)",title="t = 0.15s",titlesize=14,labelsize=12,cmap="gray",fignum=2)
                subplot(3,3,3) ; SeisPlot(u[140,:,:],pclip=100,dx=1000*dx,dy=1000*dz,xticks=[0,500,1000,1500,2000],yticks=[0,300,600,900,1200],xlabel="x (m)",ylabel="z (m)",title="t = 0.28s",titlesize=14,labelsize=12,cmap="gray",fignum=2)
                subplot(3,3,4) ; SeisPlot(u[205,:,:],pclip=100,dx=1000*dx,dy=1000*dz,xticks=[0,500,1000,1500,2000],yticks=[0,300,600,900,1200],xlabel="x (m)",ylabel="z (m)",title="t = 0.41s",titlesize=14,labelsize=12,cmap="gray",fignum=2)
                subplot(3,3,5) ; SeisPlot(u[270,:,:],pclip=100,dx=1000*dx,dy=1000*dz,xticks=[0,500,1000,1500,2000],yticks=[0,300,600,900,1200],xlabel="x (m)",ylabel="z (m)",title="t = 0.54s",titlesize=14,labelsize=12,cmap="gray",fignum=2)
                subplot(3,3,6) ; SeisPlot(u[335,:,:],pclip=100,dx=1000*dx,dy=1000*dz,xticks=[0,500,1000,1500,2000],yticks=[0,300,600,900,1200],xlabel="x (m)",ylabel="z (m)",title="t = 0.67s",titlesize=14,labelsize=12,cmap="gray",fignum=2)
                subplot(3,3,7) ; SeisPlot(u[400,:,:],pclip=100,dx=1000*dx,dy=1000*dz,xticks=[0,500,1000,1500,2000],yticks=[0,300,600,900,1200],xlabel="x (m)",ylabel="z (m)",title="t = 0.70s",titlesize=14,labelsize=12,cmap="gray",fignum=2)
                subplot(3,3,8) ; SeisPlot(u[465,:,:],pclip=100,dx=1000*dx,dy=1000*dz,xticks=[0,500,1000,1500,2000],yticks=[0,300,600,900,1200],xlabel="x (m)",ylabel="z (m)",title="t = 0.83s",titlesize=14,labelsize=12,cmap="gray",fignum=2)
                subplot(3,3,9) ; SeisPlot(u[530,:,:],pclip=100,dx=1000*dx,dy=1000*dz,xticks=[0,500,1000,1500,2000],yticks=[0,300,600,900,1200],xlabel="x (m)",ylabel="z (m)",title="t = 0.96s",titlesize=14,labelsize=12,cmap="gray",fignum=2)

    SeisPlot(d,pclip=98,dx=1,dy=dt,xlabel="Trace",ylabel="t (s)",cmap="gray",fignum=3)

end

main()
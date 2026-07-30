import csv, sys, numpy as np
csvp, win, dt = sys.argv[1], int(sys.argv[2]), float(sys.argv[3])
lines=open(csvp).readlines()
h=next(i for i,l in enumerate(lines) if l.startswith('window,'))
rows=list(csv.DictReader(lines[h:]))
w={}
for r in rows: w.setdefault(int(r['window']),[]).append(r)
side=int(round(len(w)**0.5)); stride=win//2; centre=160
tgt=[(r,c) for r in range(side) for c in range(side)
     if stride*r<=centre<stride*r+win and stride*c<=centre<stride*c+win]
best=None
for (r,c) in tgt:
    i=r*side+c; v=w[i]
    vel=np.array([float(x['vel_los_ms']) for x in v])
    dpx=np.array([float(x['disp_az_px']) for x in v])
    t=np.arange(len(vel))
    y=vel-np.polyval(np.polyfit(t,vel,1),t)
    P=np.abs(np.fft.rfft(y))**2; P[0]=0
    f=np.fft.rfftfreq(len(y),dt); k=int(np.argmax(P))
    prom=P[k]/P[1:].mean()
    rec=(f[k],prom,float(v[0]['quality']),dpx.max()-dpx.min())
    if best is None or prom>best[1]: best=rec
    print(f"    win {i:4d} q={rec[2]:.3f} p2p {rec[3]:6.2f}px  dominant {rec[0]:6.3f} Hz  prom {rec[1]:5.1f}")
hit = abs(best[0]-0.5) < 0.06
print(f"  BEST of target windows: {best[0]:.3f} Hz  prom {best[1]:.1f}  -> {'RECOVERED 0.5 Hz' if hit else 'MISSED'}")

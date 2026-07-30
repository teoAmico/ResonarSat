import csv, sys, numpy as np
def ev(p, dt, label):
    lines=open(p).readlines()
    h=next(i for i,l in enumerate(lines) if l.startswith('window,'))
    rows=list(csv.DictReader(lines[h:]))
    w={}
    for r in rows: w.setdefault(int(r['window']),[]).append(r)
    n_looks=len(next(iter(w.values()))); df=1.0/(n_looks*dt)
    best=(None,0.0); qs=[]
    for i,v in w.items():
        vel=np.array([float(x['vel_los_ms']) for x in v]); qs.append(float(v[0]['quality']))
        t=np.arange(len(vel)); y=vel-np.polyval(np.polyfit(t,vel,1),t)
        P=np.abs(np.fft.rfft(y))**2; P[0]=0
        f=np.fft.rfftfreq(len(y),dt); k=int(np.argmax(P))
        if abs(f[k]-0.5) < 2.0*df:
            prom=P[k]/P[1:].mean()
            if prom>best[1]: best=(i,prom,f[k])
    ok = best[0] is not None
    print(f"{label:52s} q_med {np.median(qs):.3f} q_max {max(qs):.3f}  "
          f"{'RECOVERED at '+format(best[2],'.3f')+' Hz prom '+format(best[1],'.1f')+f' (win {best[0]})' if ok else 'MISSED'}")
import sys
for p,dt,lab in eval(sys.argv[1]): ev(p,dt,lab)

from pathlib import Path
import csv, json, math, random, re, statistics
from openpyxl import load_workbook
from PIL import Image, ImageDraw, ImageFont

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'analysis'; OUT.mkdir(exist_ok=True)
BANDS={'LTE 800':(791e6,862e6,41.25),'900 MHz':(880e6,960e6,41.25),'LTE 1800':(1710e6,1880e6,58.3),'2100 MHz':(1920e6,2170e6,61.0),'LTE 2600':(2500e6,2690e6,61.0),'5G NR 3600':(3400e6,3800e6,61.0)}
DATE_MAP={'19_prill':('KT','19-prill-2026'),'20_prill':('KT','20-prill-2026'),'24_prill':('KT','24-prill-2026'),'26_prill':('KT','26-prill-2026'),'2_maj':('MKA','2-maj-2026'),'3_maj':('MKA','3-maj-2026'),'4_maj':('MKA','4-maj-2026'),'5_maj':('MKA','5-maj-2026')}

def identify(path):
    s=path.stem.lower().replace(' ','_')
    for key,val in DATE_MAP.items():
        if key in s:return val
    raise ValueError(path)
def point_number(s):
    m=re.search(r'P_?0*(\d+)',str(s),re.I); return int(m.group(1)) if m else None
def quantile(values,q):
    v=sorted(values); pos=(len(v)-1)*q; lo=int(pos); hi=min(lo+1,len(v)-1); f=pos-lo
    return v[lo]*(1-f)+v[hi]*f

def read_coords():
    coords={}
    for p in ROOT.rglob('Flete_Matje*.xlsx'):
        site='MKA' if 'MKA' in p.name else 'KT'; wb=load_workbook(p,read_only=True,data_only=True)
        for ws in wb.worksheets:
            rows_for_day=[]
            headers=[str(x.value).strip() if x.value is not None else '' for x in next(ws.iter_rows(min_row=1,max_row=1))]
            lat_i=next((i for i,h in enumerate(headers) if h=="Latitude'"),None)
            lon_i=next((i for i,h in enumerate(headers) if h=="Longitude'"),None)
            for row in ws.iter_rows(min_row=2,values_only=True):
                if not row or not row[0]:continue
                pn=point_number(row[0]); lat=row[lat_i] if lat_i is not None else None; lon=row[lon_i] if lon_i is not None else None
                if isinstance(lat,(int,float)) and isinstance(lon,(int,float)):
                    if not (41.30<float(lat)<41.36):
                        lat_dms_i=next((i for i,h in enumerate(headers) if h=='Latitude°'),None)
                        mm=re.search(r"(\d+)°(\d+)'([0-9.]+)",str(row[lat_dms_i] if lat_dms_i is not None else ''))
                        if mm:lat=float(mm.group(1))+float(mm.group(2))/60+float(mm.group(3))/3600
                    rows_for_day.append((float(lat),float(lon),str(row[2]),str(row[3])))
            coords[(site,ws.title)]=rows_for_day
        wb.close()
    return coords

coords=read_coords(); records=[]
for p in sorted(ROOT.rglob('Fusha_E_*.xlsx')):
    site,date=identify(p); wb=load_workbook(p,read_only=True,data_only=True)
    valid_sheets=[ws for ws in wb.worksheets if point_number(ws.title) is not None]
    for ordinal,ws in enumerate(valid_sheets):
        pn=point_number(ws.title)
        if pn is None: continue
        freq=[]; ea=[]; em=[]
        for row in ws.iter_rows(min_row=2,values_only=True):
            if isinstance(row[0],(int,float)):
                freq.append(float(row[0])); ea.append(float(row[1])); em.append(float(row[2]))
        day_coords=coords.get((site,date),[]);c=day_coords[ordinal] if ordinal<len(day_coords) else None
        rec={'site':site,'date':date,'point':pn,'sheet':ws.title,'source_file':p.name,'lat':c[0] if c else None,'lon':c[1] if c else None,'start_time':c[2] if c else '','end_time':c[3] if c else ''}
        q=0.0
        for band,(lo,hi,ref) in BANDS.items():
            idx=[i for i,f in enumerate(freq) if lo<=f<=hi]
            eact=math.sqrt(sum(ea[i]**2 for i in idx)); emax=math.sqrt(sum(em[i]**2 for i in idx))
            rec[band]=eact; rec[band+' EmaxRSS']=emax; rec[band+' ratio']=emax/eact if eact else None; q+=(eact/ref)**2
        rec['Qth_session']=q; rec['Etotal']=math.sqrt(sum(rec[b]**2 for b in BANDS)); records.append(rec)
    wb.close()
records.sort(key=lambda r:(r['site'],r['date'],r['point']))
missing=[(r['site'],r['date'],r['point'],r['sheet']) for r in records if r['lat'] is None or r['lon'] is None]
if missing:
    print('MISSING COORDS',missing)
with (OUT/'session_band_values.csv').open('w',newline='',encoding='utf-8-sig') as f:
    w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)

def xy_metres(rows):
    lon0=statistics.mean(r['lon'] for r in rows);lat0=statistics.mean(r['lat'] for r in rows);cl=math.cos(math.radians(lat0))
    return [((r['lon']-lon0)*111320*cl,(r['lat']-lat0)*110540) for r in rows]
def knn_weights(xy,k=6):
    W=[]
    for i,p in enumerate(xy):
        ds=sorted(((math.dist(p,q),j) for j,q in enumerate(xy) if j!=i))[:k];W.append([j for _,j in ds])
    return W
def moran(x,W):
    m=statistics.mean(x);z=[v-m for v in x];den=sum(v*v for v in z);num=sum(z[i]*z[j]/len(W[i]) for i in range(len(x)) for j in W[i]);return num/den
def moran_test(x,W,seed=42,nperm=999):
    obs=moran(x,W);rng=random.Random(seed);ext=0
    for _ in range(nperm):
        y=x[:];rng.shuffle(y)
        if abs(moran(y,W))>=abs(obs):ext+=1
    return obs,(ext+1)/(nperm+1)
def bh_adjust(pvals):
    n=len(pvals);order=sorted(range(n),key=lambda i:pvals[i]);adj=[1.0]*n;running=1.0
    for rank_idx in range(n-1,-1,-1):
        i=order[rank_idx];rank=rank_idx+1;running=min(running,pvals[i]*n/rank);adj[i]=running
    return adj
def local_moran(x,W,seed=42,nperm=999):
    n=len(x);m=statistics.mean(x);z=[v-m for v in x];m2=sum(v*v for v in z)/n
    def one(values,i):return values[i]*sum(values[j] for j in W[i])/len(W[i])/m2
    obs=[one(z,i) for i in range(n)];ext=[0]*n;rng=random.Random(seed)
    for _ in range(nperm):
        y=z[:];rng.shuffle(y)
        for i in range(n):
            if abs(one(y,i))>=abs(obs[i]):ext[i]+=1
    raw=[(e+1)/(nperm+1) for e in ext];adj=bh_adjust(raw);labels=[]
    for i in range(n):
        lag=sum(z[j] for j in W[i])/len(W[i])
        if adj[i]>=.05:labels.append('NS')
        elif z[i]>=0 and lag>=0:labels.append('HH')
        elif z[i]<0 and lag<0:labels.append('LL')
        elif z[i]>=0:labels.append('HL')
        else:labels.append('LH')
    return obs,raw,adj,labels
def rankdata(vals):
    order=sorted(range(len(vals)),key=lambda i:vals[i]);r=[0.0]*len(vals);i=0
    while i<len(order):
        j=i+1
        while j<len(order) and vals[order[j]]==vals[order[i]]:j+=1
        avg=(i+j-1)/2+1
        for k in range(i,j):r[order[k]]=avg
        i=j
    return r
def corr(a,b):
    ma=statistics.mean(a);mb=statistics.mean(b);da=[x-ma for x in a];db=[x-mb for x in b];den=math.sqrt(sum(x*x for x in da)*sum(x*x for x in db));return sum(x*y for x,y in zip(da,db))/den if den else 0.0
def mantel_distance_decay(x,xy,seed=42,nperm=999):
    pairs=[(i,j) for i in range(len(x)) for j in range(i+1,len(x))];dist=[math.dist(xy[i],xy[j]) for i,j in pairs];rd=rankdata(dist)
    def stat(vals):return corr(rd,rankdata([abs(vals[i]-vals[j]) for i,j in pairs]))
    obs=stat(x);rng=random.Random(seed);ext=0
    for _ in range(nperm):
        y=x[:];rng.shuffle(y)
        if abs(stat(y))>=abs(obs):ext+=1
    return obs,(ext+1)/(nperm+1)
def sign_test(a,b):
    wins=sum(x<y for x,y in zip(a,b));loss=sum(x>y for x,y in zip(a,b));n=wins+loss;k=min(wins,loss)
    p=min(1.0,2*sum(math.comb(n,i) for i in range(k+1))/(2**n)) if n else 1.0
    return wins,loss,p

summary={'n_total':len(records),'coord_missing':sum(r['lat'] is None or r['lon'] is None for r in records),'sites':{}}
for site in ('KT','MKA'):
    rows=[r for r in records if r['site']==site];xy=xy_metres(rows);W=knn_weights(xy)
    ss={'n':len(rows),'bands':{},'qth_max':max(r['Qth_session'] for r in rows),'qth_p95':quantile([r['Qth_session'] for r in rows],.95),'qth_mean':statistics.mean(r['Qth_session'] for r in rows)}
    for band in BANDS:
        x=[r[band] for r in rows];I,p=moran_test(x,W);nn=[];idw=[]
        for i,pnt in enumerate(xy):
            ds=sorted((math.dist(pnt,q),j) for j,q in enumerate(xy) if j!=i);nn.append(x[ds[0][1]])
            weights=[(1/max(d,1e-9)**2,j) for d,j in ds];idw.append(sum(w*x[j] for w,j in weights)/sum(w for w,j in weights))
        rmse_nn=math.sqrt(statistics.mean((a-b)**2 for a,b in zip(x,nn)));rmse_idw=math.sqrt(statistics.mean((a-b)**2 for a,b in zip(x,idw)))
        sensitivity={}
        for k in (4,6,8,10):
            Ik,pk=moran_test(x,knn_weights(xy,k),seed=42+k);sensitivity[str(k)]={'I':Ik,'p':pk}
        li,lp,lq,ll=local_moran(x,W,seed=84);counts={key:ll.count(key) for key in ('HH','LL','HL','LH','NS')}
        decay_r,decay_p=mantel_distance_decay(x,xy,seed=126)
        err_nn=[abs(a-b) for a,b in zip(x,nn)];err_idw=[abs(a-b) for a,b in zip(x,idw)];wins,losses,sign_p=sign_test(err_idw,err_nn)
        ss['bands'][band]={'mean':statistics.mean(x),'sd':statistics.stdev(x),'median':statistics.median(x),'p25':quantile(x,.25),'p75':quantile(x,.75),'p95':quantile(x,.95),'min':min(x),'max':max(x),'moran_I':I,'moran_p':p,'moran_sensitivity':sensitivity,'local_counts':counts,'local_labels':ll,'distance_decay_rho':decay_r,'distance_decay_p':decay_p,'rmse_nn':rmse_nn,'rmse_idw':rmse_idw,'rmse_reduction_pct':100*(rmse_nn-rmse_idw)/rmse_nn,'idw_abs_error_wins':wins,'nn_abs_error_wins':losses,'sign_p':sign_p}
    summary['sites'][site]=ss
maxrow=max(records,key=lambda r:r['Qth_session']);summary['qth_global_max']=maxrow['Qth_session'];summary['qth_global_max_row']={k:maxrow[k] for k in ('site','date','point','Qth_session')}
(OUT/'summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False),encoding='utf-8')

def hull(points):
    pts=sorted(set(points))
    def cross(o,a,b):return (a[0]-o[0])*(b[1]-o[1])-(a[1]-o[1])*(b[0]-o[0])
    lo=[]
    for p in pts:
        while len(lo)>=2 and cross(lo[-2],lo[-1],p)<=0:lo.pop()
        lo.append(p)
    up=[]
    for p in reversed(pts):
        while len(up)>=2 and cross(up[-2],up[-1],p)<=0:up.pop()
        up.append(p)
    return lo[:-1]+up[:-1]
def clip(poly,a,b,c):
    out=[]
    for P,Q in zip(poly,poly[1:]+poly[:1]):
        fp=a*P[0]+b*P[1]-c;fq=a*Q[0]+b*Q[1]-c;Pin=fp<=1e-9;Qin=fq<=1e-9
        if Pin:out.append(P)
        if Pin!=Qin:
            t=fp/(fp-fq);out.append((P[0]+t*(Q[0]-P[0]),P[1]+t*(Q[1]-P[1])))
    return out
def cells(points):
    base=hull(points);res=[]
    for i,p in enumerate(points):
        poly=base[:]
        for j,q in enumerate(points):
            if i==j:continue
            a=2*(q[0]-p[0]);b=2*(q[1]-p[1]);c=q[0]*q[0]+q[1]*q[1]-p[0]*p[0]-p[1]*p[1];poly=clip(poly,a,b,c)
            if not poly:break
        res.append(poly)
    return res,base

palette=['#440154','#3b528b','#21918c','#5ec962','#fde725']
try: font=ImageFont.truetype('arial.ttf',22);small=ImageFont.truetype('arial.ttf',16);tiny=ImageFont.truetype('arial.ttf',13)
except: font=small=tiny=ImageFont.load_default()
all_by_band={b:[r[b] for r in records] for b in BANDS}
for site in ('KT','MKA'):
    rows=[r for r in records if r['site']==site];xy=xy_metres(rows);polys,bound=cells(xy)
    img=Image.new('RGB',(2400,1600),'white');d=ImageDraw.Draw(img);d.text((60,20),f'Voronoi partition of band-integrated RF-EMF measurements at {site}',fill='black',font=font)
    for idx,band in enumerate(BANDS):
        col=idx%3;row=idx//3;x0=60+col*780;y0=90+row*730;W=600;H=600;xs=[p[0] for p in xy];ys=[p[1] for p in xy];xmin,xmax=min(xs),max(xs);ymin,ymax=min(ys),max(ys)
        def conv(p):return (x0+(p[0]-xmin)/(xmax-xmin)*W,y0+H-(p[1]-ymin)/(ymax-ymin)*H)
        vals=all_by_band[band];cuts=[quantile(vals,q) for q in (.2,.4,.6,.8)]
        for poly,r in zip(polys,rows):
            k=sum(r[band]>c for c in cuts);d.polygon([conv(p) for p in poly],fill=palette[k],outline='white')
        d.line([conv(p) for p in bound+[bound[0]]],fill='black',width=3)
        for p in xy:
            xp,yp=conv(p);d.ellipse((xp-3,yp-3,xp+3,yp+3),fill='black')
        d.text((x0,y0-28),band,fill='black',font=small);d.text((x0+W-25,y0+15),'N',fill='black',font=small);d.line((x0+W-16,y0+70,x0+W-16,y0+35),fill='black',width=3);d.polygon([(x0+W-16,y0+28),(x0+W-24,y0+42),(x0+W-8,y0+42)],fill='black')
        px50=50/(xmax-xmin)*W;d.line((x0+20,y0+H-20,x0+20+px50,y0+H-20),fill='black',width=4);d.text((x0+20,y0+H-45),'50 m',fill='black',font=tiny)
        for k,color in enumerate(palette):
            xx=x0+W+15;yy=y0+40+k*32;d.rectangle((xx,yy,xx+22,yy+22),fill=color);label=(f'≤ {cuts[0]:.2f}' if k==0 else f'> {cuts[-1]:.2f}' if k==4 else f'{cuts[k-1]:.2f}-{cuts[k]:.2f}');d.text((xx+28,yy),label,fill='black',font=tiny)
        d.text((x0+W+15,y0+10),'Eact RSS (V/m)',fill='black',font=tiny)
    img.save(OUT/f'voronoi_{site}.png',dpi=(300,300))

cluster_palette={'HH':'#d73027','LL':'#4575b4','HL':'#fdae61','LH':'#74add1','NS':'#d9d9d9'}
for site in ('KT','MKA'):
    rows=[r for r in records if r['site']==site];xy=xy_metres(rows);polys,bound=cells(xy)
    img=Image.new('RGB',(2400,1600),'white');d=ImageDraw.Draw(img);d.text((60,20),f'Local Moran cluster maps at {site} (BH-FDR q < 0.05)',fill='black',font=font)
    for idx,band in enumerate(BANDS):
        col=idx%3;row=idx//3;x0=60+col*780;y0=90+row*730;W=600;H=600;xs=[p[0] for p in xy];ys=[p[1] for p in xy];xmin,xmax=min(xs),max(xs);ymin,ymax=min(ys),max(ys)
        def conv(p):return (x0+(p[0]-xmin)/(xmax-xmin)*W,y0+H-(p[1]-ymin)/(ymax-ymin)*H)
        labels=summary['sites'][site]['bands'][band]['local_labels']
        for poly,label in zip(polys,labels):d.polygon([conv(p) for p in poly],fill=cluster_palette[label],outline='white')
        d.line([conv(p) for p in bound+[bound[0]]],fill='black',width=3)
        for p in xy:
            xp,yp=conv(p);d.ellipse((xp-3,yp-3,xp+3,yp+3),fill='black')
        d.text((x0,y0-28),band,fill='black',font=small);d.text((x0+W-25,y0+15),'N',fill='black',font=small);d.line((x0+W-16,y0+70,x0+W-16,y0+35),fill='black',width=3);d.polygon([(x0+W-16,y0+28),(x0+W-24,y0+42),(x0+W-8,y0+42)],fill='black')
        for k,label in enumerate(('HH','LL','HL','LH','NS')):
            xx=x0+W+15;yy=y0+40+k*32;d.rectangle((xx,yy,xx+22,yy+22),fill=cluster_palette[label]);d.text((xx+28,yy),label,fill='black',font=tiny)
    img.save(OUT/f'local_moran_{site}.png',dpi=(300,300))

print(json.dumps({'records':len(records),'by_site':{s:sum(r['site']==s for r in records) for s in ('KT','MKA')},'qth_max':summary['qth_global_max'],'coord_missing':summary['coord_missing']},indent=2))

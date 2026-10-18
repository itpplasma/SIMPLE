"""Retrieve bibliographic metadata from Crossref; read-only network research."""
import requests,json,concurrent.futures
from pathlib import Path
known={'burby':'10.1063/1.5119801','pauli':'10.1016/j.cpc.2021.107981','dvi':'10.1063/1.5022277','piecewise':'10.1051/m2an:2008006','simple':'10.1088/1361-6587/ae10dd','carybrizard':'10.1103/RevModPhys.81.693','brizardhahm':'10.1103/RevModPhys.79.421','ricketson':'10.1016/j.jcp.2020.109639','filtered':'10.1007/s00211-020-01105-3','large':'10.1007/s00211-022-01298-9','ua3d':'10.1137/19M127402X','wangzhao':'10.1137/20M1340101','mrcm':'10.1007/s00211-013-0602-0','dimits':'10.1063/1.3683000','speiser':'10.1029/JZ070i017p04219','buchner':'10.1029/JA094iA09p11821','caryescande':'10.1103/PhysRevA.34.4256','tennyson':'10.1103/PhysRevLett.56.2117','gni':'10.1007/3-540-30666-8','markl':'10.3217/efp2p-0x485'}
queries={'stroboscopic':'Numerical stroboscopic averaging for ODEs and DAEs Calvo Chartier Murua Sanz-Serna','higheraverage':'Higher-Order Averaging Formal Series Numerical Integration I B-series Chartier Murua Sanz-Serna','driftlorentz':'An implicit drift Lorentz mover for arbitrarily magnetized species Cohen Friedman Vay','envelope':'Envelope following method dynamical systems Petzold','twoexp':'Two-scale exponential integrators with uniform accuracy for three-dimensional charged-particle dynamics under strong magnetic field','ua2d':'Uniformly accurate methods for Vlasov equations with non-homogeneous strong magnetic field Chartier Crouseilles Lemou Mehats Zhao'}
def get(k,v,isquery=False):
 try:
  url='https://api.crossref.org/works';r=requests.get(url if isquery else url+'/'+v,params={'query.bibliographic':v,'rows':2} if isquery else {},timeout=40);r.raise_for_status();m=r.json()['message'];return k,{'source':r.url,'records':m['items'] if isquery else [m]}
 except Exception as e:return k,{'error':str(e)}
with concurrent.futures.ThreadPoolExecutor(max_workers=6) as ex:
 fut=[ex.submit(get,k,v) for k,v in known.items()]+[ex.submit(get,k,v,True) for k,v in queries.items()]
 out=dict(f.result() for f in concurrent.futures.as_completed(fut))
Path(__file__).with_name('citation_metadata.json').write_text(json.dumps(out,indent=2))
for k,d in sorted(out.items()):
 if 'error' in d:print(k,d);continue
 for m in d['records'][:1]:print(k,m.get('title'),[(a.get('given'),a.get('family')) for a in m.get('author',[])],m.get('published'),m.get('container-title'),m.get('volume'),m.get('page'),m.get('article-number'),m.get('DOI'))

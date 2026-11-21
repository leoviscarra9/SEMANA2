
const dot=(a,b)=>a.reduce((s,v,i)=>s+v*b[i],0),nrm=a=>Math.sqrt(dot(a,a)),mv=(A,x)=>A.map(r=>dot(r,x)),sub=(a,b)=>a.map((v,i)=>v-b[i]);
function luSolve(A0,b0){const n=b0.length,U=A0.map(r=>r.slice()),b=b0.slice(),L=Array.from({length:n},(_,i)=>Array(n).fill(0).map((_,j)=>i==j?1:0));
for(let k=0;k<n-1;k++){let m=k;for(let i=k+1;i<n;i++)if(Math.abs(U[i][k])>Math.abs(U[m][k]))m=i;
if(m!=k){[U[k],U[m]]=[U[m],U[k]];[b[k],b[m]]=[b[m],b[k]];for(let j=0;j<k;j++)[L[k][j],L[m][j]]=[L[m][j],L[k][j]];}
for(let i=k+1;i<n;i++){const l=U[i][k]/U[k][k];L[i][k]=l;for(let j=k;j<n;j++)U[i][j]-=l*U[k][j];}}
const y=Array(n).fill(0),x=Array(n).fill(0);for(let i=0;i<n;i++){let s=b[i];for(let j=0;j<i;j++)s-=L[i][j]*y[j];y[i]=s;}
for(let i=n-1;i>=0;i--){let s=y[i];for(let j=i+1;j<n;j++)s-=U[i][j]*x[j];x[i]=s/U[i][i];}return{x,L,U};}
function iter(kind,A,b,tol,mx){const n=b.length;let x=Array(n).fill(0),h=[],nb=nrm(b),r=b.slice(),p=b.slice(),rs=dot(r,r);
for(let k=1;k<=mx;k++){
 if(kind=='J'){x=x.map((xi,i)=>(b[i]-dot(A[i],x)+A[i][i]*xi)/A[i][i]);}
 else if(kind=='G'){for(let i=0;i<n;i++){let s=b[i];for(let j=0;j<n;j++)if(j!=i)s-=A[i][j]*x[j];x[i]=s/A[i][i];}}
 else{const Ap=mv(A,p),al=rs/dot(p,Ap);x=x.map((v,i)=>v+al*p[i]);r=r.map((v,i)=>v-al*Ap[i]);const rn=dot(r,r);p=r.map((v,i)=>v+rn/rs*p[i]);rs=rn;}
 const res=nrm(sub(b,mv(A,x)))/nb;h.push(res);
 if(!isFinite(res))return{x,it:k,ok:false,hist:h,div:true};
 if(res<tol)return{x,it:k,ok:true,hist:h};}
return{x,it:mx,ok:false,hist:h};}
function cond2(A){const n=A.length;let v=Array(n).fill(0).map((_,i)=>1+0.1*i),lm=0;
for(let k=0;k<600;k++){const w=mv(A,v),nn=nrm(w);v=w.map(t=>t/nn);lm=dot(v,mv(A,v));}
let u=Array(n).fill(1),lmin=0;for(let k=0;k<600;k++){const w=luSolve(A,u).x,nn=nrm(w);u=w.map(t=>t/nn);lmin=1/dot(u,luSolve(A,u).x);}
return lm/lmin;}
const SYS={
 NORMAL:{A:[[10,-1,2,0],[-1,11,-1,3],[2,-1,10,-1],[0,3,-1,8]],xs:[1,2,-1,1],nota:'SPD, diagonal dominante, n=4'},
 ESTRES:{A:(()=>{const n=8;return Array.from({length:n},(_,i)=>Array.from({length:n},(_,j)=>i==j?2:Math.abs(i-j)==1?-1:0))})(),xs:[1,2,3,4,5,6,7,8],nota:'SPD tridiagonal 2,−1; n=8'},
 MAL:{A:Array.from({length:5},(_,i)=>Array.from({length:5},(_,j)=>1/(i+j+1))),xs:[1,1,1,1,1],nota:'Matriz de Hilbert 5×5 (SPD, κ elevado)'}};
function runLin(esc,tol,mx){const S=SYS[esc],A=S.A,xs=S.xs,b=mv(A,xs),nx=nrm(xs),out=[];
const ev=(name,x,it,ok,hist,div)=>({name,it,ok,div:!!div,res:nrm(sub(b,mv(A,x)))/nrm(b),eabs:nrm(sub(x,xs)),erel:nrm(sub(x,xs))/nx,hist:hist||[],x});
let r=iter('J',A,b,tol,mx);out.push(ev('Jacobi',r.x,r.it,r.ok,r.hist,r.div));
r=iter('G',A,b,tol,mx);out.push(ev('Gauss–Seidel',r.x,r.it,r.ok,r.hist,r.div));
const l=luSolve(A,b);out.push(ev('LU (pivoteo parcial)',l.x,1,true,[]));
r=iter('C',A,b,tol,mx);out.push(ev('Gradiente conjugado',r.x,r.it,r.ok,r.hist,r.div));
const db=b.map((_,i)=>1e-8*nrm(b)/Math.sqrt(b.length)*(i%2?1:-1)),xp=luSolve(A,b.map((v,i)=>v+db[i])).x;
return{A,b,xs,kappa:cond2(A),rows:out,sens:{relb:nrm(db)/nrm(b),relx:nrm(sub(xp,xs))/nx}};}
function model(p){
 if(p.esc=='MAL'){const c=p.c,T=p.ts;return{v0:c*T*T,acc:t=>2*c*(T-t),v:t=>c*(T-t)**2,vp:t=>-2*c*(T-t),x:t=>c*(T**3-(T-t)**3)/3,td:T,dbl:true};}
 const{v0,a,k}=p;return{v0,acc:(t,v)=>a+k*v,v:t=>(v0+a/k)*Math.exp(-k*t)-a/k,vp:t=>-(v0+a/k)*k*Math.exp(-k*t),x:t=>(v0+a/k)*(1-Math.exp(-k*t))/k-a*t/k,td:Math.log(1+v0*k/a)/k};}
function rk(M,order,h){const n=Math.max(1,Math.round(M.td/h)),hh=M.td/n,f=(t,y)=>[y[1],-M.acc(t,y[1])];let t=0,y=[0,M.v0];const T=[0],X=[0],V=[M.v0];
for(let i=0;i<n;i++){let k1=f(t,y);
 if(order==2){const k2=f(t+hh,[y[0]+hh*k1[0],y[1]+hh*k1[1]]);y=[y[0]+hh/2*(k1[0]+k2[0]),y[1]+hh/2*(k1[1]+k2[1])];}
 else{const a=(c,s)=>[y[0]+s*c[0],y[1]+s*c[1]],k2=f(t+hh/2,a(k1,hh/2)),k3=f(t+hh/2,a(k2,hh/2)),k4=f(t+hh,a(k3,hh));y=[0,1].map(j=>y[j]+hh/6*(k1[j]+2*k2[j]+2*k3[j]+k4[j]));}
 t+=hh;T.push(t);X.push(y[0]);V.push(y[1]);}
const ev=Math.max(...T.map((tt,i)=>Math.abs(V[i]-M.v(tt)))),ex=Math.max(...T.map((tt,i)=>Math.abs(X[i]-M.x(tt))));
return{order,n,h:hh,T,X,V,errV:ev,errX:ex,errXT:Math.abs(X[n]-M.x(M.td)),xT:X[n]};}
function quad(M){const ex=M.x(M.td),rows=[],tr=N=>{const h=M.td/N;let s=0.5*(M.v(0)+M.v(M.td));for(let i=1;i<N;i++)s+=M.v(i*h);return h*s;},
sp=N=>{const h=M.td/N;let s=M.v(0)+M.v(M.td);for(let i=1;i<N;i++)s+=(i%2?4:2)*M.v(i*h);return h/3*s;};
const R=[];for(let i=0;i<5;i++){R[i]=[tr(2**i)];for(let j=1;j<=i;j++)R[i][j]=(4**j*R[i][j-1]-R[i-1][j-1])/(4**j-1);}
for(let i=2;i<=5;i++){const N=2**i;rows.push({N,trap:tr(N),eT:Math.abs(tr(N)-ex),simp:sp(N),eS:Math.abs(sp(N)-ex)});}
return{ex,rows,romb:R.map((r,i)=>({N:2**i,val:r[i],err:Math.abs(r[i]-ex)}))};}
function roots(M,tol=1e-10){const g=M.v,dg=M.vp,td=M.td,out=[];
{let lo=0,hi=2*td;if(g(lo)*g(hi)>0)out.push({m:'Bisección',na:'sin cambio de signo en [0, 2t*]: no aplicable'});else{let it=0;while(hi-lo>tol&&it<200){const m=(lo+hi)/2;it++;if(g(lo)*g(m)<=0)hi=m;else lo=m;}const t=(lo+hi)/2;out.push({m:'Bisección',it,t,err:Math.abs(t-td)});}}
{let t=M.dbl?0.6*td:0,it=0;for(;it<500;){if(g(t)===0)break;const tn=t-g(t)/dg(t);it++;if(!isFinite(tn))break;if(Math.abs(tn-t)<tol){t=tn;break;}t=tn;}out.push({m:'Newton–Raphson',it,t,err:Math.abs(t-td)});}
if(M.dbl){let t=0.6*td,it=0;for(;it<50;){if(g(t)===0)break;const tn=t-2*g(t)/dg(t);it++;if(!isFinite(tn))break;if(Math.abs(tn-t)<tol){t=tn;break;}t=tn;}out.push({m:'Newton modificado (m=2)',it,t,err:Math.abs(t-td)});}
{let a=M.dbl?0.6*td:0,b=M.dbl?0.8*td:1.5*td,it=0;for(;it<500;){const d=g(b)-g(a);if(d==0)break;const c=b-g(b)*(b-a)/d;it++;a=b;b=c;if(Math.abs(b-a)<tol)break;}out.push({m:'Secante',it,t:b,err:Math.abs(b-td)});}
return out;}
function sens1(p){const M=model(p),r={};
if(p.esc=='MAL'){const e=1e-6;r.eps=e;r.shift=Math.sqrt(e/p.c);r.amp=r.shift/e;}
else{const M2=model({...p,a:p.a*1.01}),M3=model({...p,v0:p.v0*1.01});r.da=(M2.td-M.td)/M.td/0.01;r.dv=(M3.td-M.td)/M.td/0.01;}
return r;}
function sim1(p){const M=model(p),r2=rk(M,2,p.h),r4=rk(M,4,p.h),q=quad(M),hs=[4,2,1,.5,.25].map(f=>p.h*f);
const conv=hs.map(h=>({h,e2:rk(M,2,h).errXT,e4:rk(M,4,h).errXT}));
const xd=M.x(M.td);return{M,td:M.td,dfren:xd,dtot:xd+M.v0*p.tr,rk2:r2,rk4:r4,q,roots:roots(M),conv,sens:sens1(p)};}
const $=s=>document.querySelector(s),F=(x,d=3)=>x===undefined?'—':!isFinite(x)?String(x):x===0?'0':(Math.abs(x)<1e-3||Math.abs(x)>=1e5)?x.toExponential(d):x.toFixed(6);
const PRE={NORMAL:{esc:'NORMAL',v0:5,a:.5,k:.05,h:.5,tr:1,c:.05,ts:10},ESTRES:{esc:'ESTRES',v0:6,a:.35,k:.08,h:2,tr:2,c:.05,ts:10},MAL:{esc:'MAL',v0:5,a:.5,k:.05,h:.5,tr:1,c:.05,ts:10}};
const CL=['#df765b','#087e83','#537b91','#b58b43'];
const LBL={NORMAL:'NORMAL',ESTRES:'ESTRÉS',MAL:'MAL CONDICIONADO'};
const MODS=[
{id:'m1',t:'Detención',ic:'🛑',st:'Operativo (S1)',q:'¿Cómo varían tiempo y distancia de detención ante cambios de velocidad, reacción y frenado?',mod:"dx/dt = v,  dv/dt = −(a + k·v).  Referencia: v(t) = (v₀ + a/k)e^(−kt) − a/k;  t_det = ln(1 + v₀k/a)/k.",met:'RK2 (Heun), RK4, trapecio, Simpson 1/3, Romberg, bisección, Newton, secante'},
{id:'m2',t:'Deterioro',ic:'🔧',st:'Sistemas lineales operativos (S1)',q:'¿Cuándo alcanzaría un indicador hipotético de deterioro el umbral de inspección?',mod:"D(t) = D₀ + αt + βt²;  resolver D(t) = D_umbral. Calibración: sistema SPD A·x = b.",met:'Jacobi, Gauss–Seidel, LU, gradiente conjugado; raíces (S2)'},
{id:'m3',t:'Oscilaciones',ic:'🌬️',st:'Modelo definido; métodos en S2',q:'¿Cómo cambia la oscilación ante fuerza externa, amortiguamiento y parámetros?',mod:"x′ = v,  v′ = (F(t) − c·v − k·x)/m.  Referencia: F=0, subamortiguado, x(t)=e^(−ζω₀t)[C₁cos ω_d t + C₂ sen ω_d t].",met:'RK2, RK4, trapecio, Simpson, Romberg'},
{id:'m4',t:'Sensores',ic:'📡',st:'Modelo definido; métodos en S2',q:'¿Cómo reconstruir registros faltantes y medir la distorsión por ruido?',mod:"Datos sintéticos y(t)=sen t o e^(−t); se ocultan puntos y se mide el error real.",met:'Lagrange, Newton (diferencias divididas), spline cúbico, raíces'}];
const ESC3=[['NORMAL','Parámetros nominales, paso adecuado'],['ESTRES','Paso grueso, mayor velocidad y reacción'],['MAL','Raíz doble (v=c(t*−t)²) / matriz de Hilbert']];
const H={get(){try{return JSON.parse(localStorage.getItem('simtele_hist')||'[]')}catch(e){return[]}},set(a){try{localStorage.setItem('simtele_hist',JSON.stringify(a))}catch(e){}}};
function save(caso,p,res){const a=H.get();a.unshift({id:Date.now(),fecha:new Date().toLocaleString('es-BO'),caso,esc:p.esc,p,res});H.set(a.slice(0,50));}
function chart(lines,o={}){const W=360,Hh=220,m={l:46,r:10,t:10,b:34},pts=lines.flatMap(l=>l.pts.filter(q=>isFinite(q[0])&&isFinite(q[1])&&(!o.ly||q[1]>0)&&(!o.lx||q[0]>0)));
if(!pts.length)return'<p class=note>Sin puntos que graficar.</p>';const tx=o.lx?Math.log10:x=>x,ty=o.ly?Math.log10:x=>x;
let xs=pts.map(q=>tx(q[0])),ys=pts.map(q=>ty(q[1])),x0=Math.min(...xs),x1=Math.max(...xs),y0=Math.min(...ys),y1=Math.max(...ys);if(o.ly){y0=Math.floor(y0);y1=Math.ceil(y1)}if(x1==x0)x1++;if(y1==y0)y1++;
const X=v=>m.l+(tx(v)-x0)/(x1-x0)*(W-m.l-m.r),Y=v=>Hh-m.b-(ty(v)-y0)/(y1-y0)*(Hh-m.t-m.b);let s=`<svg viewBox="0 0 ${W} ${Hh}" role=img aria-label="${o.title||''}"><rect x=${m.l} y=${m.t} width=${W-m.l-m.r} height=${Hh-m.t-m.b} class=pf />`;
for(let i=0;i<=4;i++){const yv=y0+(y1-y0)*i/4,xv=x0+(x1-x0)*i/4,yl=o.ly?'1e'+Math.round(yv):(+yv.toPrecision(2)),xl=o.lx?(10**xv).toPrecision(2):(+xv.toPrecision(2));
s+=`<text x=${m.l-4} y=${Hh-m.b-(i/4)*(Hh-m.t-m.b)+3} text-anchor=end class=tk>${yl}</text><text x=${m.l+(i/4)*(W-m.l-m.r)} y=${Hh-m.b+12} text-anchor=middle class=tk>${xl}</text>`;}
lines.forEach((l,i)=>{const p=l.pts.filter(q=>isFinite(q[0])&&isFinite(q[1])&&(!o.ly||q[1]>0)&&(!o.lx||q[0]>0));s+=`<polyline fill=none class="s${i}" ${l.dash?'stroke-dasharray="4 3"':''} points="${p.map(q=>X(q[0]).toFixed(1)+','+Y(q[1]).toFixed(1)).join(' ')}"/>`;if(l.dots)p.forEach(q=>s+=`<circle cx=${X(q[0]).toFixed(1)} cy=${Y(q[1]).toFixed(1)} r=2.5 class="f${i}"/>`);});
s+=`<text x=${(W+m.l)/2} y=${Hh-3} text-anchor=middle class=ax>${o.xl||''}</text><text transform="rotate(-90)" x=${-Hh/2} y=11 text-anchor=middle class=ax>${o.yl||''}</text></svg><div class=leg>${lines.map((l,i)=>`<span><i class="k${i}"></i>${l.name}</span>`).join('')}</div>`;return s;}
function tbl(h,rows){return`<div class=tw><table><thead><tr>${h.map(x=>`<th>${x}</th>`).join('')}</tr></thead><tbody>${rows.map(r=>`<tr>${r.map(c=>`<td>${c}</td>`).join('')}</tr>`).join('')}</tbody></table></div>`}
function home(){return`<div class=hero><h1>Simula. Compara. Verifica.</h1><p>Cuatro casos de estudio, tres escenarios y un solo objetivo: entender cuánto podemos confiar en un número. Plataforma educativa inspirada en un teleférico andino; datos sintéticos.</p><div class=chips><span>4 módulos</span><span>3 escenarios</span><span>13 métodos</span><span>errores y condicionamiento</span></div><svg class=hsvg viewBox="0 0 700 150" preserveAspectRatio="none"><defs><linearGradient id=mg x1=0 y1=0 x2=0 y2=1><stop offset=0 stop-color="#6d9da2"/><stop offset=1 stop-color="#456e7b"/></linearGradient><linearGradient id=mg2 x1=0 y1=0 x2=0 y2=1><stop offset=0 stop-color="#315b69" stop-opacity=".9"/><stop offset=1 stop-color="#1d3d4d"/></linearGradient></defs><g><circle class=star cx=60 cy=20 r=1.5 /><circle class=star cx=180 cy=45 r=1.2 /><circle class=star cx=330 cy=15 r=1.6 /><circle class=star cx=470 cy=40 r=1.2 /><circle class=star cx=600 cy=18 r=1.5 /><circle class=star cx=660 cy=55 r=1.1 /></g><path d="M0 150V95L70 60L130 98L215 40L300 105L390 55L470 100L560 35L640 90L700 65V150Z" fill="url(#mg)"/><path d="M0 150V125L90 92L170 125L260 85L350 128L450 90L540 125L620 95L700 120V150Z" fill="url(#mg2)"/><line x1=0 y1=38 x2=700 y2=-30 stroke="#d8eeee" stroke-opacity=".6" stroke-width="1.4"/><line x1=0 y1=48 x2=700 y2=-20 stroke="#d8eeee" stroke-opacity=".3" stroke-width="1"/><g class=gond><line x1=0 y1=0 x2=0 y2=12 stroke="#d8eeee" stroke-width="1.2" transform="translate(40,-14)"/><rect x=26 y=-2 width=28 height=18 rx=5 fill="#df765b"/><rect x=31 y=3 width=8 height=7 rx=1.5 fill="#fff" opacity=".85"/><rect x=42 y=3 width=8 height=7 rx=1.5 fill="#fff" opacity=".85"/></g><g class=gond2 transform="translate(0,10)"><rect x=26 y=-2 width=28 height=18 rx=5 fill="#087e83"/><rect x=31 y=3 width=8 height=7 rx=1.5 fill="#fff" opacity=".85"/><rect x=42 y=3 width=8 height=7 rx=1.5 fill="#fff" opacity=".85"/></g></svg></div><div class=grid>${MODS.map((m,i)=>`<a class=card style="--a:${CL[i]}" href="#${m.id}"><span class=ic>${m.ic}</span><h3>${m.t}</h3><p>${m.q}</p><span class=badge>${m.st}</span></a>`).join('')}</div><div class=grid2><a class=card style="--a:#0b6e99" href="#hist"><h3>Historial de experimentos</h3><p>Guardados en este navegador; abrir y relanzar.</p></a><a class=card style="--a:#e4572e" href="#err"><h3>Teoría de errores (transversal)</h3><p>Definiciones usadas en todos los módulos.</p></a></div>`}
function modInfo(m){return`<h1 class=mh style="--a:${CL[MODS.indexOf(m)]}">${m.ic} Módulo: ${m.t}</h1><p class=lead>${m.q}</p><div class=box><b>Modelo</b><br>${m.mod}<br><b>Métodos</b>: ${m.met}</div><h3>Escenarios definidos</h3>${tbl(['Escenario','Definición'],ESC3.map(e=>[LBL[e[0]],e[1]]))}`}
let cur1='NORMAL';
function m1(){const p=PRE[cur1],isM=cur1=='MAL';return modInfo(MODS[0])+`<div class=ctl><label>Escenario<select id=esc>${ESC3.map(e=>`<option value=${e[0]} ${e[0]==cur1?'selected':''}>${LBL[e[0]]}</option>`).join('')}</select></label>`+
(isM?[['c','c (m/s³)'],['ts','t* (s)']]:[['v0','v₀ (m/s)'],['a','a (m/s²)'],['k','k (1/s)']]).concat([['h','h (s)'],['tr','t_reacc (s)']]).map(f=>`<label>${f[1]}<input id=${f[0]} type=number step=any value=${p[f[0]]}></label>`).join('')+`<button id=run>Ejecutar</button></div><div id=out></div>`}
function readP(){const p={esc:$('#esc').value};for(const k of['v0','a','k','h','tr','c','ts']){const e=$('#'+k);p[k]=e?parseFloat(e.value):PRE[p.esc][k];}return p;}
function valid(p){const e=[];if(p.esc=='MAL'){if(!(p.c>0))e.push('c > 0');if(!(p.ts>0))e.push('t* > 0')}else{if(!(p.v0>0))e.push('v₀ > 0');if(!(p.a>0))e.push('a > 0');if(!(p.k>0))e.push('k > 0')}if(!(p.h>0))e.push('h > 0');if(!(p.tr>=0))e.push('t_reacc ≥ 0');if(p.h>(p.esc=='MAL'?p.ts:1e9))e.push('h ≤ t*');return e;}
function show1(p,s){const M=s.M,T=Array.from({length:101},(_,i)=>i/100*s.td);
$('#out').innerHTML=`<h3>Resultados</h3>${tbl(['Magnitud','Valor'],[['t_det (analítico)',F(s.td)+' s'],['Distancia de frenado (analítica)',F(s.dfren)+' m'],['Distancia total (+ v₀·t_reacc)',F(s.dtot)+' m'],['Pasos RK (h ajustado)',s.rk4.n+' (h = '+F(s.rk4.h)+' s)']])}
<h3>EDO: RK2 vs RK4 frente a la solución analítica</h3>${tbl(['Método','Máx |error v|','Máx |error x|','|error x(t_det)|','% rel. x'],[s.rk2,s.rk4].map(r=>['RK'+r.order+(r.order==2?' (Heun)':''),F(r.errV),F(r.errX),F(r.errXT),F(100*r.errXT/s.dfren)]))}
<div class=two><div>${chart([{name:'Exacta',pts:T.map(t=>[t,M.v(t)])},{name:'RK2',pts:s.rk2.T.map((t,i)=>[t,s.rk2.V[i]]),dash:1,dots:1},{name:'RK4',pts:s.rk4.T.map((t,i)=>[t,s.rk4.V[i]]),dots:1}],{xl:'t (s)',yl:'v (m/s)',title:'Velocidad'})}</div><div>${chart([{name:'RK2',pts:s.conv.map(c=>[c.h,c.e2]),dots:1},{name:'RK4',pts:s.conv.map(c=>[c.h,c.e4]),dots:1}],{lx:1,ly:1,xl:'h (s)',yl:'|error x(t_det)| (m)',title:'Convergencia'})}</div></div>
<h3>Integración de v(t) sobre [0, t_det] (= distancia)</h3>${tbl(['N','Trapecio','Error','Simpson 1/3','Error'],s.q.rows.map(r=>[r.N,F(r.trap),F(r.eT),F(r.simp),F(r.eS)]))}${tbl(['Romberg R(i,i)','N','Valor','Error'],s.q.romb.map((r,i)=>[i,r.N,F(r.val),F(r.err)]))}
<h3>Raíces de v(t)=0</h3>${tbl(['Método','Iteraciones','t estimado','|error|'],s.roots.map(r=>r.na?[r.m,'—','—',r.na]:[r.m,r.it,F(r.t,8),F(r.err)]))}
<h3>Sensibilidad</h3><p>${s.sens.eps?`Perturbar v en ε = ${s.sens.eps}: las raíces se desplazan ${F(s.sens.shift)} s (√(ε/c)); amplificación ≈ ${F(s.sens.amp,0)}×. Es mal condicionamiento de una <b>raíz doble</b>, no un fallo de convergencia aleatorio.`:`Un +1 % en <i>a</i> cambia t_det en ${F(100*s.sens.da)} % y un +1 % en v₀ lo cambia en ${F(100*s.sens.dv)} % (número de condición relativo ≈ ${F(Math.abs(s.sens.da))} y ${F(s.sens.dv)}).`}</p>
<h3>Interpretación</h3><div class=two><div class=box><b>Conclusión numérica</b><br>${s.rk4.errXT<s.rk2.errXT?'RK4 reduce el error de RK2 en un factor ≈ '+F(s.rk2.errXT/Math.max(s.rk4.errXT,1e-300),2)+' con el mismo h.':'RK2 y RK4 alcanzan precisión similar.'} ${p.esc=='MAL'?'Bisección no aplica (sin cambio de signo); Newton converge linealmente (factor ½); la corrección por multiplicidad m=2 restaura la convergencia rápida.':'Simpson/Romberg superan al trapecio en orden de convergencia; Newton y secante convergen en pocas iteraciones.'}</div><div class=box><b>Conclusión del problema</b><br>Con estos parámetros sintéticos la cabina se detiene en ≈ ${F(s.td,4)} s y recorre ≈ ${F(s.dtot,4)} m incluyendo la reacción. No se infiere seguridad real: el modelo ignora pendiente, carga y fallos de equipo.</div></div><p class=note>Experimento guardado en el historial.</p>`;}
function run1(){const p=readP(),e=valid(p);if(e.length){$('#out').innerHTML='<p class=err>Parámetros inválidos: '+e.join(', ')+'</p>';return}const s=sim1(p);show1(p,s);save('Detención',p,{td:s.td,dfren:s.dfren,dtot:s.dtot,e2:s.rk2.errXT,e4:s.rk4.errXT});}
let cur2='NORMAL';
function m2(){return modInfo(MODS[1])+`<h3>Sistemas lineales A·x = b (A SPD, x* conocida)</h3><div class=ctl><label>Escenario<select id=esc2>${ESC3.map(e=>`<option value=${e[0]} ${e[0]==cur2?'selected':''}>${LBL[e[0]]}</option>`).join('')}</select></label><label>Tolerancia<input id=tol type=number step=any value=1e-10></label><label>Máx. iter.<input id=mx type=number value=20000></label><button id=run2>Ejecutar</button></div><div id=out></div>`}
function run2(){const esc=$('#esc2').value,tol=parseFloat($('#tol').value),mx=parseInt($('#mx').value);if(!(tol>0&&mx>0)){$('#out').innerHTML='<p class=err>Tolerancia y máx. iter. deben ser > 0</p>';return}const r=runLin(esc,tol,mx),S=SYS[esc];
$('#out').innerHTML=`<p><b>${LBL[esc]}</b>: ${S.nota}. κ₂(A) ≈ ${F(r.kappa)}. x* = [${S.xs.join(', ')}]. Criterio: ‖b−Ax‖/‖b‖ &lt; tol.</p>${tbl(['Método','Iter.','Convergió','Residuo rel.','Error rel.','Error abs.'],r.rows.map(x=>[x.name,x.it,x.div?'diverge':x.ok?'sí':'no (máx. iter.)',F(x.res),F(x.erel),F(x.eabs)]))}
${chart(r.rows.filter(x=>x.hist.length).map((x,i)=>({name:x.name,pts:x.hist.map((v,k)=>[k+1,v])})),{ly:1,xl:'Iteración',yl:'Residuo relativo',title:'Residuo'})}
<p>Perturbación relativa de b de ${F(r.sens.relb)} produce un cambio relativo de x de ${F(r.sens.relx)} (amplificación ${F(r.sens.relx/r.sens.relb,2)}× ≤ κ ≈ ${F(r.kappa)}).</p>
<div class=two><div class=box><b>Conclusión numérica</b><br>${esc=='MAL'?'Jacobi diverge (la matriz no es diagonalmente dominante); Gauss–Seidel converge pero lento. LU y CG dan el error esperado ≈ κ·ε_máq. El mal condicionamiento ≠ fallo de convergencia.':'CG converge en ≤ n iteraciones (aritmética exacta); Jacobi y GS convergen por dominancia o SPD, GS en ≈ la mitad de iteraciones.'}</div><div class=box><b>Conclusión del problema</b><br>Estos sistemas calibran coeficientes del modelo de deterioro. Un κ alto implica que pequeños errores de medición se amplifican en los coeficientes estimados.</div></div>`;
const p={esc,tol,mx};save('Deterioro (sist. lineales)',p,{kappa:r.kappa,rows:r.rows.map(x=>({n:x.name,it:x.it,res:x.res,er:x.erel}))});}
function errs(){return`<h1>Teoría de errores (transversal)</h1>${tbl(['Concepto','Definición','Uso'],[['Error absoluto','E_a = |x − x̃|','Todos los módulos'],['Error relativo','E_r = |x − x̃| / |x| (si x ≠ 0)','Se omite si x=0'],['Error porcentual','100·E_r','Tablas de resumen'],['Residual','‖b − A x̃‖/‖b‖ ; |f(x̃)|','Sistemas lineales y raíces'],['Sensibilidad','κ(A) = ‖A‖‖A⁻¹‖ ;  Δx/x ≤ κ Δb/b','Módulos 1 y 2 (S1)']])}<p class=note>Residuo pequeño no implica error pequeño si κ es grande.</p>`}
function hist(){const a=H.get();return`<h1>Historial de experimentos</h1>${a.length?tbl(['Fecha','Caso','Escenario','Acciones'],a.map(r=>[r.fecha,r.caso,LBL[r.esc],`<button class=sm data-o=${r.id}>Abrir</button> <button class=sm data-r=${r.id}>Relanzar</button> <button class=sm data-d=${r.id}>Borrar</button>`])):'<p>Aún no hay experimentos. Ejecute un caso.</p>'}<div id=out></div><p class=note>Persistencia en el navegador (localStorage). El esquema SQL relacional está en el repositorio (ver informe).</p>`}
function route(){const h=(location.hash||'#home').slice(1);document.querySelectorAll('nav a').forEach(a=>a.classList.toggle('on',a.getAttribute('href')=='#'+h));const v=$('#view');
v.innerHTML=h=='m1'?m1():h=='m2'?m2():h=='m3'?modInfo(MODS[2]):h=='m4'?modInfo(MODS[3]):h=='hist'?hist():h=='err'?errs():home();
if(h=='m1'){$('#esc').onchange=e=>{cur1=e.target.value;v.innerHTML=m1();route2()};$('#run').onclick=run1;}
if(h=='m2'){$('#esc2').onchange=e=>cur2=e.target.value;$('#run2').onclick=run2;}
if(h=='hist')v.onclick=hclick;else v.onclick=null;}
function route2(){$('#esc').onchange=e=>{cur1=e.target.value;$('#view').innerHTML=m1();route2()};$('#run').onclick=run1;}
function hclick(e){const t=e.target,a=H.get();if(t.dataset.d){H.set(a.filter(r=>r.id!=t.dataset.d));route();return}
const r=a.find(r=>r.id==(t.dataset.o||t.dataset.r));if(!r)return;
if(t.dataset.o)$('#out').innerHTML=`<h3>Experimento guardado</h3><pre>${JSON.stringify({caso:r.caso,escenario:r.esc,parametros:r.p,resultados:r.res},null,1)}</pre>`;
if(t.dataset.r){let ok,msg;if(r.caso=='Detención'){const s=sim1(r.p);ok=Math.abs(s.td-r.res.td)<1e-9&&Math.abs(s.dtot-r.res.dtot)<1e-9&&Math.abs(s.rk4.errXT-r.res.e4)<1e-9;msg=`t_det: ${F(r.res.td,8)} → ${F(s.td,8)}; d_total: ${F(r.res.dtot,8)} → ${F(s.dtot,8)}`}else{const l=runLin(r.p.esc,r.p.tol,r.p.mx);ok=Math.abs(l.kappa-r.res.kappa)/r.res.kappa<1e-9&&l.rows.every((x,i)=>x.it==r.res.rows[i].it);msg=`κ: ${F(r.res.kappa)} → ${F(l.kappa)}; iteraciones coinciden: ${ok}`}
$('#out').innerHTML=`<h3>Relanzamiento</h3><p class=${ok?'ok':'err'}>${ok?'✔ Coincide':'✘ No coincide'} dentro de la tolerancia 1e-9. ${msg}</p>`;}}
window.addEventListener('hashchange',route);route();

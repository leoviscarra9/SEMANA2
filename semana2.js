/* SIMTELE-RISK · Semana 2 (extiende script.js): interpolación, raíces, M3 oscilaciones, M4 sensores, comparación y base de datos IndexedDB */
MODS[0].st='Operativo (S1–S2)';MODS[1].st='Lineales + raíces (S1–S2)';MODS[1].met+=', bisección, Newton, secante';MODS[2].st='Operativo (S2)';MODS[3].st='Operativo (S2)';
MODS[3].met='Lagrange, Newton (diferencias divididas), spline cúbico natural, bisección, Newton, secante';
const mi=(m,defs)=>modInfo(m).replace(/<h3>Escenarios definidos[\s\S]*$/,'<h3>Escenarios definidos</h3>'+tbl(['Escenario','Definición'],defs.map(e=>[LBL[e[0]],e[1]])));
const er=x=>x.na?-1:x.err;
/* ---------- Raíces (genérico) ---------- */
function rt(m,g,dg,a,b,t0,tol,ex,q){const h=[];let it=0,t;
if(m=='B'){let lo=a,hi=b;if(g(lo)*g(hi)>0)return{m:'Bisección',na:'sin cambio de signo en [a, b]'};while(hi-lo>tol&&it<200){const c=(lo+hi)/2;it++;if(g(lo)*g(c)<=0)hi=c;else lo=c;h.push(Math.abs((lo+hi)/2-ex))}t=(lo+hi)/2;return{m:'Bisección',it,t,err:Math.abs(t-ex),h}}
if(m=='N'){t=t0;q=q||1;for(;it<500;){const d=dg(t);if(!d||!isFinite(d))break;const tn=t-q*g(t)/d;if(!isFinite(tn))break;it++;h.push(Math.abs(tn-ex));const ok=Math.abs(tn-t)<tol;t=tn;if(ok)break}return{m:q>1?'Newton modificado (m=2)':'Newton–Raphson',it,t,err:Math.abs(t-ex),h}}
let x0=t0,x1=t0+.1*(b-a);for(;it<500;){const d=g(x1)-g(x0);if(!d||!isFinite(d))break;const c=x1-g(x1)*(x1-x0)/d;if(!isFinite(c))break;it++;x0=x1;x1=c;h.push(Math.abs(c-ex));if(Math.abs(x1-x0)<tol)break}return{m:'Secante',it,t:x1,err:Math.abs(x1-ex),h}}
const rtab=r=>tbl(['Método','Iteraciones','t estimado','|error|'],r.map(x=>x.na?[x.m,'—','—',x.na]:[x.m,x.it,F(x.t,8),F(x.err)]));
const rchart=r=>chart(r.filter(x=>x.h&&x.h.length).map(x=>({name:x.m,pts:x.h.map((e,k)=>[k+1,e])})),{ly:1,xl:'Iteración',yl:'|error| de la raíz',title:'Convergencia de raíces'});
/* ---------- M2: deterioro D(t)=D0+αt+βt² ---------- */
const DET={NORMAL:{D0:1,al:.8,be:.02,Du:60,d:'D(t) creciente, raíz simple'},ESTRES:{D0:1,al:.001,be:.02,Du:1e5,d:'D′(0)≈0 y umbral muy alto: Newton desde t₀=0 salta lejos'},MAL:{D0:0,al:2,be:-.05,Du:20,d:'Umbral = máximo de D(t): raíz doble en t*=20'}};
function det(esc){const{D0,al,be,Du}=DET[esc],g=t=>D0+al*t+be*t*t-Du,dg=t=>al+2*be*t,ex=(-al+Math.sqrt(Math.max(0,al*al-4*be*(D0-Du))))/(2*be),M=esc=='MAL',t0=M?.6*ex:0,hi=2*ex,tol=1e-10;
const r=[rt('B',g,dg,0,hi,0,tol,ex),rt('N',g,dg,0,0,t0,tol,ex)];if(M)r.push(rt('N',g,dg,0,0,t0,tol,ex,2));r.push(rt('S',g,dg,0,hi,t0,tol,ex));return{D0,al,be,Du,ex,r}}
const flatDet=esc=>{const r=det(esc).r;return{eB:er(r[0]),eN:er(r[1]),eS:er(r[r.length-1]),itN:r[1].it}};
function ext2(){const esc=$('#esc2').value,d=det(esc),M=esc=='MAL',T=Array.from({length:101},(_,i)=>i/100*1.3*d.ex);
$('#out').insertAdjacentHTML('beforeend',`<h3>Raíces: instante en que D(t) alcanza el umbral (S2)</h3><p>D(t) = ${d.D0} + ${d.al}·t + (${d.be})·t² ; D<sub>u</sub> = ${d.Du} (${DET[esc].d}). Raíz exacta (fórmula cuadrática): t* = ${F(d.ex,8)}.</p>${rtab(d.r)}
<div class=two><div>${chart([{name:'D(t)',pts:T.map(t=>[t,d.D0+d.al*t+d.be*t*t])},{name:'Umbral',pts:[[0,d.Du],[T[100],d.Du]],dash:1}],{xl:'t',yl:'D(t)',title:'Deterioro'})}</div><div>${rchart(d.r)}</div></div>
<div class=box><b>Conclusión numérica</b><br>${M?'Bisección no aplica (D−D<sub>u</sub> no cambia de signo). Newton converge linealmente (error ×½ por iteración); con m=2 recupera convergencia cuadrática. La raíz estimada queda limitada a ≈√ε_máq por mal condicionamiento.':esc=='ESTRES'?'Con D′(0)≈0 el primer paso de Newton cae en t≈10⁸ y luego converge monótonamente desde la derecha: '+d.r[1].it+' iteraciones frente a pocas del caso NORMAL.':'Bisección converge linealmente (≈'+d.r[0].it+' it.); Newton y secante cuadrática/superlinealmente.'}</div>`);
save('Deterioro (raíces)',{esc},flatDet(esc))}
/* ---------- M3: oscilador x'' = −(c x' + k x)/m ---------- */
const OSC={NORMAL:{m:1,c:.4,k:4,x0:1,v0:0,T:10,h:.1},ESTRES:{m:1,c:.1,k:4,x0:1,v0:0,T:30,h:.45},MAL:{m:1,c:200,k:4,x0:1,v0:0,T:2,h:.02}};
const D3=[['NORMAL','Subamortiguado (ζ=0,1), h=0,1 s, T=10 s'],['ESTRES','Amortiguamiento bajo, T=30 s y paso grueso h=0,45 s'],['MAL','Sistema rígido (c=200, λ≈−200): h=0,02 s excede la región de estabilidad']];
function oscEx(p){const{m,c,k,x0,v0}=p,d=c*c-4*m*k;
if(d<0){const a=-c/(2*m),w=Math.sqrt(-d)/(2*m),C=(v0-a*x0)/w;return{x:t=>Math.exp(a*t)*(x0*Math.cos(w*t)+C*Math.sin(w*t)),v:t=>Math.exp(a*t)*((a*x0+w*C)*Math.cos(w*t)+(a*C-w*x0)*Math.sin(w*t))}}
const s=Math.sqrt(d),r1=(-c+s)/(2*m),r2=(-c-s)/(2*m),A=(v0-r2*x0)/(r1-r2),B=x0-A;return{x:t=>A*Math.exp(r1*t)+B*Math.exp(r2*t),v:t=>A*r1*Math.exp(r1*t)+B*r2*Math.exp(r2*t)}}
function rkO(p,E,o,h){const n=Math.max(1,Math.round(p.T/h)),hh=p.T/n,f=y=>[y[1],(-p.c*y[1]-p.k*y[0])/p.m];let y=[p.x0,p.v0],ex=0,ev=0;const T=[0],X=[y[0]];
for(let i=0;i<n;i++){const k1=f(y);if(o==2){const k2=f([y[0]+hh*k1[0],y[1]+hh*k1[1]]);y=[y[0]+hh/2*(k1[0]+k2[0]),y[1]+hh/2*(k1[1]+k2[1])]}else{const a=(c,s)=>[y[0]+s*c[0],y[1]+s*c[1]],k2=f(a(k1,hh/2)),k3=f(a(k2,hh/2)),k4=f(a(k3,hh));y=[0,1].map(j=>y[j]+hh/6*(k1[j]+2*k2[j]+2*k3[j]+k4[j]))}
const t=(i+1)*hh;T.push(t);X.push(y[0]);ex=Math.max(ex,Math.abs(y[0]-E.x(t)));ev=Math.max(ev,Math.abs(y[1]-E.v(t)))}
return{o,n,h:hh,T,X,ex,ev,eT:Math.abs(y[0]-E.x(p.T)),xT:y[0]}}
function quadG(f,T,ex){const tr=N=>{const h=T/N;let s=.5*(f(0)+f(T));for(let i=1;i<N;i++)s+=f(i*h);return h*s},sp=N=>{const h=T/N;let s=f(0)+f(T);for(let i=1;i<N;i++)s+=(i%2?4:2)*f(i*h);return h/3*s};
const R=[];for(let i=0;i<5;i++){R[i]=[tr(2**i)];for(let j=1;j<=i;j++)R[i][j]=(4**j*R[i][j-1]-R[i-1][j-1])/(4**j-1)}
return{ex,rows:[2,3,4,5].map(i=>{const N=2**i;return{N,trap:tr(N),eT:Math.abs(tr(N)-ex),simp:sp(N),eS:Math.abs(sp(N)-ex)}}),romb:R.map((r,i)=>({N:2**i,val:r[i],err:Math.abs(r[i]-ex)}))}}
function sim3(p){const E=oscEx(p),en=(x,v)=>.5*p.m*v*v+.5*p.k*x*x,xT=E.x(p.T),en0=en(p.x0,p.v0),W=en0-en(xT,E.v(p.T)),r2=rkO(p,E,2,p.h),r4=rkO(p,E,4,p.h),
conv=[4,2,1,.5,.25].map(f=>p.h*f).map(h=>({h,e2:rkO(p,E,2,h).eT,e4:rkO(p,E,4,h).eT})),q=quadG(t=>p.c*E.v(t)**2,p.T,W);return{E,xT,en0,W,r2,r4,conv,q}}
const flat3=s=>({xT:s.xT,W:s.W,e2:s.r2.eT,e4:s.r4.eT,eS:s.q.rows[3].eS,eR:s.q.romb[4].err});
let cur3='NORMAL';
function m3(){const p=OSC[cur3];return mi(MODS[2],D3)+`<div class=ctl><label>Escenario<select id=e3>${D3.map(e=>`<option value=${e[0]} ${e[0]==cur3?'selected':''}>${LBL[e[0]]}</option>`).join('')}</select></label>`+[['m','m (kg)'],['c','c (N·s/m)'],['k','k (N/m)'],['x0','x₀ (m)'],['v0','v₀ (m/s)'],['T','T (s)'],['h','h (s)']].map(f=>`<label>${f[1]}<input id=i_${f[0]} type=number step=any value=${p[f[0]]}></label>`).join('')+`<button id=r3>Ejecutar</button></div><div id=out></div>`}
function run3(){const p={esc:$('#e3').value};for(const k of['m','c','k','x0','v0','T','h'])p[k]=parseFloat($('#i_'+k).value);const e=[];if(!(p.m>0))e.push('m > 0');if(!(p.k>0))e.push('k > 0');if(!(p.c>=0))e.push('c ≥ 0');if(!(p.T>0))e.push('T > 0');if(!(p.h>0&&p.h<=p.T))e.push('0 < h ≤ T');if(!isFinite(p.x0+p.v0))e.push('x₀, v₀ numéricos');if(e.length){$('#out').innerHTML='<p class=err>Parámetros inválidos: '+e.join(', ')+'</p>';return}
const s=sim3(p),T=Array.from({length:201},(_,i)=>i/200*p.T),B=10*(Math.abs(p.x0)+1),cl=v=>Math.abs(v)<B?v:NaN,rig=p.c*p.c-4*p.m*p.k>=0,bad=!isFinite(s.r4.eT)||s.r4.eT>1;
$('#out').innerHTML=`<h3>Resultados</h3>${tbl(['Magnitud','Valor'],[['Tipo',rig?'sobreamortiguado (rígido)':'subamortiguado'],['x(T) analítica',F(s.xT)+' m'],['Energía inicial E₀',F(s.en0)+' J'],['Energía disipada ∫c·v²dt = E₀−E(T)',F(s.W)+' J'],['Pasos RK (h ajustado)',s.r4.n+' (h = '+F(s.r4.h)+' s)']])}
<h3>EDO: RK2 vs RK4 frente a la solución analítica</h3>${tbl(['Método','Máx |error x|','Máx |error v|','|error x(T)|','% rel. x(T)'],[s.r2,s.r4].map(r=>['RK'+r.o+(r.o==2?' (Heun)':''),F(r.ex),F(r.ev),F(r.eT),F(100*r.eT/Math.abs(s.xT))]))}
<div class=two><div>${chart([{name:'Exacta',pts:T.map(t=>[t,s.E.x(t)])},{name:'RK2',pts:s.r2.T.map((t,i)=>[t,cl(s.r2.X[i])]),dash:1},{name:'RK4',pts:s.r4.T.map((t,i)=>[t,cl(s.r4.X[i])]),dots:1}],{xl:'t (s)',yl:'x (m)',title:'Posición'})}</div><div>${chart([{name:'RK2',pts:s.conv.map(c=>[c.h,c.e2]),dots:1},{name:'RK4',pts:s.conv.map(c=>[c.h,c.e4]),dots:1}],{lx:1,ly:1,xl:'h (s)',yl:'|error x(T)| (m)',title:'Convergencia'})}</div></div>
<h3>Integración de la potencia disipada c·v(t)² sobre [0, T]</h3>${tbl(['N','Trapecio','Error','Simpson 1/3','Error'],s.q.rows.map(r=>[r.N,F(r.trap),F(r.eT),F(r.simp),F(r.eS)]))}${tbl(['Romberg R(i,i)','N','Valor','Error'],s.q.romb.map((r,i)=>[i,r.N,F(r.val),F(r.err)]))}
<div class=two><div class=box><b>Conclusión numérica</b><br>${bad?'RK2 y RK4 <b>divergen/son inestables</b>: h·|λ| queda fuera de la región de estabilidad explícita (RK4 ≈ 2,78). Al reducir h por debajo de ≈2,78/|λ| el error cae de golpe (ver convergencia): no es un problema de orden sino de estabilidad.':'RK4 reduce el error de RK2 en ≈'+F(s.r2.eT/Math.max(s.r4.eT,1e-300),2)+'× con el mismo h; Simpson y Romberg superan al trapecio al integrar la energía disipada.'}</div><div class=box><b>Conclusión del problema</b><br>Con parámetros sintéticos, la energía disipada es ≈ ${F(s.W,4)} J en ${p.T} s. No se infiere seguridad estructural real.</div></div><p class=note>Experimento guardado en el historial.</p>`;
save('Oscilaciones',p,flat3(s))}
/* ---------- M4: interpolación ---------- */
const lag=(X,Y,x)=>{let s=0;for(let i=0;i<X.length;i++){let L=Y[i];for(let j=0;j<X.length;j++)if(j!=i)L*=(x-X[j])/(X[i]-X[j]);s+=L}return s};
function nd(X,Y){const n=X.length,c=Y.slice();for(let j=1;j<n;j++)for(let i=n-1;i>=j;i--)c[i]=(c[i]-c[i-1])/(X[i]-X[i-j]);return x=>{let s=c[n-1];for(let i=n-2;i>=0;i--)s=s*(x-X[i])+c[i];return s}}
function spl(X,Y){const n=X.length,h=[],a=Y,c=Array(n).fill(0),b=[],d=[],l=[1],mu=[0],z=[0];for(let i=0;i<n-1;i++)h[i]=X[i+1]-X[i];
for(let i=1;i<n-1;i++){const al=3/h[i]*(a[i+1]-a[i])-3/h[i-1]*(a[i]-a[i-1]);l[i]=2*(X[i+1]-X[i-1])-h[i-1]*mu[i-1];mu[i]=h[i]/l[i];z[i]=(al-h[i-1]*z[i-1])/l[i]}
for(let j=n-2;j>=0;j--){c[j]=z[j]-mu[j]*c[j+1];b[j]=(a[j+1]-a[j])/h[j]-h[j]*(c[j+1]+2*c[j])/3;d[j]=(c[j+1]-c[j])/(3*h[j])}
return x=>{let i=0;while(i<n-2&&x>X[i+1])i++;const t=x-X[i];return a[i]+b[i]*t+c[i]*t*t+d[i]*t*t*t}}
const SG={NORMAL:{f:Math.sin,a:0,b:2*Math.PI,n:9,r:0,rt:{l:.5,a:0,b:Math.PI/2,t0:.3,ex:Math.PI/6}},ESTRES:{f:t=>Math.exp(-t)*Math.cos(3*t),a:0,b:4,n:7,r:.05,rt:{l:0,a:0,b:1,t0:.3,ex:Math.PI/6}},MAL:{f:t=>1/(1+25*t*t),a:-1,b:1,n:15,r:0,rt:{l:.5,a:0,b:1,t0:.1,ex:.2}}};
const D4=[['NORMAL','y = sen t en [0, 2π], 9 nodos, sin ruido'],['ESTRES','y = e^(−t)·cos 3t en [0, 4], 7 nodos y ruido del 5 %'],['MAL','Runge 1/(1+25t²) en [−1, 1], 15 nodos equidistantes (Runge)']];
function sim4(p){const S=SG[p.esc],f=S.f,{a,b}=S,n=p.n,X=Array.from({length:n},(_,i)=>a+(b-a)*i/(n-1)),Yt=X.map(f);let sd=12345;const rn=()=>(sd=(sd*1664525+1013904223)%4294967296)/4294967296-.5,A=Math.max(...Yt.map(Math.abs)),Y=Yt.map(y=>y+2*p.r*A*rn()),
PN=nd(X,Y),SP=spl(X,Y),PL=t=>lag(X,Y,t),G=Array.from({length:201},(_,i)=>a+(b-a)*i/200),ev=P=>{let m=0,s=0;G.forEach(t=>{const e=Math.abs(P(t)-f(t));m=Math.max(m,e);s+=e*e});return{max:m,rms:Math.sqrt(s/G.length)}},
eL=ev(PL),eN=ev(PN),eS=ev(SP),dLN=Math.max(...G.map(t=>Math.abs(PL(t)-PN(t)))),R=S.rt,g=t=>SP(t)-R.l,dg=t=>(SP(t+1e-6)-SP(t-1e-6))/2e-6,
roots=[rt('B',g,dg,R.a,R.b,0,1e-10,R.ex),rt('N',g,dg,0,0,R.t0,1e-10,R.ex),rt('S',g,dg,R.a,R.b,R.t0,1e-10,R.ex)];return{f,X,Y,G,PL,PN,SP,eL,eN,eS,dLN,roots,R}}
const flat4=s=>({eL:s.eL.max,eN:s.eN.max,eS:s.eS.max,dLN:s.dLN,rB:er(s.roots[0]),rN:er(s.roots[1]),rS:er(s.roots[2])});
let cur4='NORMAL';
function m4(){const S=SG[cur4];return mi(MODS[3],D4)+`<div class=ctl><label>Escenario<select id=e4>${D4.map(e=>`<option value=${e[0]} ${e[0]==cur4?'selected':''}>${LBL[e[0]]}</option>`).join('')}</select></label><label>Nodos n<input id=i_n type=number value=${S.n}></label><label>Ruido (fracción)<input id=i_r type=number step=any value=${S.r}></label><button id=r4>Ejecutar</button></div><div id=out></div>`}
function run4(){const p={esc:$('#e4').value,n:parseInt($('#i_n').value),r:parseFloat($('#i_r').value)},e=[];if(!(p.n>=3&&p.n<=40))e.push('3 ≤ n ≤ 40');if(!(p.r>=0&&p.r<=.5))e.push('0 ≤ ruido ≤ 0,5');if(e.length){$('#out').innerHTML='<p class=err>Parámetros inválidos: '+e.join(', ')+'</p>';return}
const s=sim4(p),G=s.G,ln=(name,P,d)=>({name,pts:G.map(t=>[t,P(t)]),dash:d});
$('#out').innerHTML=`<h3>Interpolación frente a la función verdadera (201 puntos de evaluación)</h3>${tbl(['Método','Máx |error|','Error RMS'],[['Lagrange',F(s.eL.max),F(s.eL.rms)],['Newton (dif. divididas)',F(s.eN.max),F(s.eN.rms)],['Spline cúbico natural',F(s.eS.max),F(s.eS.rms)]])}<p class=note>Máx |Lagrange − Newton| = ${F(s.dLN)} (mismo polinomio; la diferencia es error de redondeo).</p>
<div class=two><div>${chart([ln('Exacta',s.f),ln('Lagrange/Newton',s.PN,1),ln('Spline',s.SP),{name:'Datos',pts:s.X.map((x,i)=>[x,s.Y[i]]),dots:1}],{xl:'t',yl:'y',title:'Interpolación'})}</div><div>${chart([{name:'Lagrange',pts:G.map(t=>[t,Math.abs(s.PL(t)-s.f(t))])},{name:'Newton',pts:G.map(t=>[t,Math.abs(s.PN(t)-s.f(t))]),dash:1},{name:'Spline',pts:G.map(t=>[t,Math.abs(s.SP(t)-s.f(t))])}],{ly:1,xl:'t',yl:'|error|',title:'Error'})}</div></div>
<h3>Raíces sobre el spline: s(t) = ${s.R.l} (valor exacto t* = ${F(s.R.ex,8)})</h3>${rtab(s.roots)}
<div class=two><div class=box><b>Conclusión numérica</b><br>${p.esc=='MAL'?'El polinomio de grado '+(p.n-1)+' oscila en los extremos (fenómeno de Runge, error máx '+F(s.eL.max,2)+'); el spline lo evita ('+F(s.eS.max,2)+'). Más nodos equidistantes empeoran el polinomio.':p.esc=='ESTRES'?'Con ruido, interpolar exactamente copia el ruido: el error ya no baja con n. El spline es más estable que el polinomio.':'Con datos limpios, Lagrange y Newton coinciden y el spline logra el menor error en los extremos.'}</div><div class=box><b>Conclusión del problema</b><br>Reconstruir registros de sensor faltantes con un polinomio global es arriesgado con muchos nodos o con ruido; conviene spline. Los datos son sintéticos.</div></div><p class=note>Experimento guardado en el historial.</p>`;
save('Sensores',p,flat4(s))}
/* ---------- Comparación de escenarios ---------- */
function cmp(){const E=['NORMAL','ESTRES','MAL'],rows=[];E.forEach(e=>{const s=sim1(PRE[e]),l=runLin(e,1e-10,20000),d=det(e),s3=sim3({esc:e,...OSC[e]}),p4={esc:e,n:SG[e].n,r:SG[e].r},s4=sim4(p4),L=LBL[e],bad=l.rows.filter(x=>!x.ok).map(x=>x.name+(x.div?' diverge':' no converge'));
rows.push(['1 Detención',L,`RK2 ${F(s.rk2.errXT)} · RK4 ${F(s.rk4.errXT)}`,`Simpson ${F(s.q.rows[3].eS)} · Romberg ${F(s.q.romb[4].err)}`,s.roots.some(x=>x.na)?'Bisección no aplica':'ok']);
rows.push(['2 Deterioro',L,`κ=${F(l.kappa)} · it J/GS/CG: ${l.rows[0].it}/${l.rows[1].it}/${l.rows[3].it}`,`Raíz: Bis ${F(er(d.r[0]))} · Newton ${F(er(d.r[1]))} (${d.r[1].it} it)`,bad.join('; ')||(d.r[0].na?'Bisección no aplica':'ok')]);
rows.push(['3 Oscilaciones',L,`RK2 ${F(s3.r2.eT)} · RK4 ${F(s3.r4.eT)}`,`Simpson ${F(s3.q.rows[3].eS)} · Romberg ${F(s3.q.romb[4].err)}`,!isFinite(s3.r4.eT)||s3.r4.eT>1?'RK inestable':'ok']);
rows.push(['4 Sensores',L,`Lagr. ${F(s4.eL.max)} · Newton ${F(s4.eN.max)} · spline ${F(s4.eS.max)}`,`Raíz: Bis ${F(er(s4.roots[0]))} · Newton ${F(er(s4.roots[1]))} · Sec ${F(er(s4.roots[2]))}`,s4.eL.max>1?'Runge: polinomio oscila':'ok'])});
return`<h1>Comparación de escenarios</h1><p class=lead>Primeras comparaciones NORMAL / ESTRÉS / MAL CONDICIONADO en los cuatro módulos (errores máximos o en el extremo; tolerancia 1e-10).</p>${tbl(['Módulo','Escenario','Indicador 1','Indicador 2','Observación'],rows)}<p class=note>Lectura: el error crece de NORMAL a ESTRÉS por paso/ruido y en MAL por condicionamiento o estabilidad, no por fallos del código.</p>`}
/* ---------- Base de datos (IndexedDB) ---------- */
let DB=null,cache=[];
const idb=(m,fn)=>new Promise(r=>{try{const t=DB.transaction('ejecuciones',m),q=fn(t.objectStore('ejecuciones'));t.oncomplete=()=>r(q&&q.result!==undefined?q.result:true);t.onerror=()=>r(null)}catch(e){r(null)}});
const dbOpen=()=>new Promise(r=>{try{const q=indexedDB.open('simtele_db',1);q.onupgradeneeded=()=>{const s=q.result.createObjectStore('ejecuciones',{keyPath:'id'});s.createIndex('caso','caso');s.createIndex('esc','esc')};q.onsuccess=()=>{DB=q.result;r(true)};q.onerror=()=>r(false);q.onblocked=()=>r(false)}catch(e){r(false)}});
const dbAll=()=>idb('readonly',s=>s.getAll()),dbGet=id=>idb('readonly',s=>s.get(id)),dbPut=o=>idb('readwrite',s=>s.put(o)),dbDel=id=>idb('readwrite',s=>s.delete(id));
H.get=()=>cache;
H.set=a=>{const o=new Set(cache.map(r=>r.id)),n=new Set(a.map(r=>r.id));cache=a;if(!DB){try{localStorage.setItem('simtele_hist',JSON.stringify(a))}catch(e){}return}a.forEach(r=>{if(!o.has(r.id))dbPut(r)});o.forEach(id=>{if(!n.has(id))dbDel(id)})};
const REL={'Detención':p=>{const s=sim1(p);return{td:s.td,dfren:s.dfren,dtot:s.dtot,e2:s.rk2.errXT,e4:s.rk4.errXT}},'Deterioro (sist. lineales)':p=>{const l=runLin(p.esc,p.tol,p.mx);return{kappa:l.kappa}},'Deterioro (raíces)':p=>flatDet(p.esc),'Oscilaciones':p=>flat3(sim3(p)),'Sensores':p=>flat4(sim4(p))};
hclick=function(e){const t=e.target,a=H.get();if(t.dataset.d){H.set(a.filter(r=>r.id!=t.dataset.d));route();return}
const r=a.find(r=>r.id==(t.dataset.o||t.dataset.r));if(!r)return;
if(t.dataset.o){$('#out').innerHTML=`<h3>Experimento guardado</h3><pre>${JSON.stringify({caso:r.caso,escenario:r.esc,parametros:r.p,resultados:r.res},null,1)}</pre>`;return}
const f=REL[r.caso],n=f?f(r.p):{},bad=Object.keys(r.res).filter(k=>typeof r.res[k]=='number').filter(k=>!(Object.is(n[k],r.res[k])||Math.abs(n[k]-r.res[k])<=1e-9*Math.max(1,Math.abs(r.res[k])))),ok=f&&!bad.length;
$('#out').innerHTML=`<h3>Relanzamiento</h3><p class=${ok?'ok':'err'}>${ok?'✔ Coincide':'✘ No coincide'} con lo guardado (tolerancia 1e-9). ${bad.length?'Difieren: '+bad.join(', '):''}</p>`};
const q1=s=>"'"+String(s).replace(/'/g,"''")+"'",dl=(n,x,ty)=>{const a=document.createElement('a');a.href=URL.createObjectURL(new Blob([x],{type:ty}));a.download=n;a.click()};
const sqlOut=()=>'-- SIMTELE-RISK: volcado del historial (esquema en schema.sql)\n'+H.get().map(r=>`INSERT INTO ejecucion(id,fecha,modulo,escenario,parametros) VALUES(${r.id},${q1(r.fecha)},${q1(r.caso)},${q1(r.esc)},${q1(JSON.stringify(r.p))});\n`+Object.entries(r.res).filter(([k,v])=>typeof v=='number'&&isFinite(v)).map(([k,v])=>`INSERT INTO metrica(ejecucion_id,nombre,valor) VALUES(${r.id},${q1(k)},${v});`).join('\n')).join('\n');
async function selfTest(){const p={...PRE.NORMAL},s=sim1(p),rec={id:Date.now(),fecha:new Date().toLocaleString('es-BO'),caso:'Detención',esc:'NORMAL',p,res:{td:s.td,dfren:s.dfren,dtot:s.dtot,e2:s.rk2.errXT,e4:s.rk4.errXT}};
let got;if(DB){await dbPut(rec);got=await dbGet(rec.id)}else{H.set([rec,...cache]);got=JSON.parse(localStorage.getItem('simtele_hist')||'[]').find(r=>r.id==rec.id)}
cache=[rec,...cache.filter(r=>r.id!=rec.id)];const s2=got&&sim1(got.p),same=!!got&&JSON.stringify(got)===JSON.stringify(rec),re=!!got&&Math.abs(s2.td-got.res.td)<1e-12&&Math.abs(s2.dtot-got.res.dtot)<1e-12;
route();$('#out2').innerHTML=`<h3>Prueba guardar → recuperar</h3>${tbl(['Paso','Resultado'],[['1. Guardar ejecución #'+rec.id+' en '+(DB?'IndexedDB':'localStorage'),'✔'],['2. Leer de la base (no de memoria)',got?'✔':'✘'],['3. Registro recuperado idéntico al guardado',same?'✔':'✘'],['4. Relanzar con los parámetros recuperados reproduce t_det y d_total',re?'✔':'✘']])}<p class=${same&&re?'ok':'err'}>${same&&re?'Prueba superada.':'Prueba fallida.'}</p>`}
function ext3(){$('#view').insertAdjacentHTML('beforeend',`<h3>Base de datos</h3><div class=ctl><button id=tst>Probar guardar y recuperar</button><button class=sm id=sqlb>Mostrar SQL</button><button class=sm id=jsb>Descargar JSON</button><button class=sm id=sqd>Descargar .sql</button></div><p class=note>Motor: ${DB?'IndexedDB «simtele_db», almacén «ejecuciones»':'localStorage (IndexedDB no disponible)'} · ${H.get().length} ejecuciones guardadas.</p><div id=out2></div>`);
$('#tst').onclick=selfTest;$('#sqlb').onclick=()=>$('#out2').innerHTML='<pre>'+sqlOut().replace(/</g,'&lt;')+'</pre>';$('#jsb').onclick=()=>dl('simtele_historial.json',JSON.stringify(H.get(),null,1),'application/json');$('#sqd').onclick=()=>dl('simtele_historial.sql',sqlOut(),'text/plain')}
/* ---------- Ruteo ---------- */
const route0=route;window.removeEventListener('hashchange',route0);
route=function(){const h=(location.hash||'#home').slice(1);if(!['m3','m4','cmp'].includes(h)){route0();if(h=='m2')$('#run2').onclick=()=>{run2();ext2()};if(h=='hist')ext3();return}
document.querySelectorAll('nav a').forEach(a=>a.classList.toggle('on',a.getAttribute('href')=='#'+h));const v=$('#view');v.onclick=null;
if(h=='cmp'){v.innerHTML=cmp();return}
const mount=()=>{v.innerHTML=h=='m3'?m3():m4();const s=$(h=='m3'?'#e3':'#e4');s.onchange=e=>{if(h=='m3')cur3=e.target.value;else cur4=e.target.value;mount()};$(h=='m3'?'#r3':'#r4').onclick=h=='m3'?run3:run4};mount()};
window.addEventListener('hashchange',route);
dbOpen().then(async ok=>{if(ok){cache=((await dbAll())||[]).sort((x,y)=>y.id-x.id);if(!cache.length){let o=[];try{o=JSON.parse(localStorage.getItem('simtele_hist')||'[]')}catch(e){}for(const r of o)await dbPut(r);cache=o}}else{try{cache=JSON.parse(localStorage.getItem('simtele_hist')||'[]')}catch(e){cache=[]}}route()});

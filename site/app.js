/* Pareto Atlas — deployed site. Vanilla JS; talks to /api/{vote,leaderboard}. */
(function(){
"use strict";
var D = window.APPDATA;
if(!D){document.body.innerHTML='<p style="padding:40px;font-family:sans-serif">Data failed to load.</p>';return;}
var PI=Math.PI, HP=PI/2, SQRT2=Math.SQRT2, SQRT3=Math.sqrt(3), XT=D.xterms, YT=D.yterms, NX=D.nx;

/* ---------- projection engine ---------- */
function familyFn(a,b){return function(lat,lon){var u=lon/PI,v=lat/HP,x=0,y=0,k,i,j;
  for(k=0;k<XT.length;k++){i=XT[k][0];j=XT[k][1];x+=a[k]*Math.pow(u,2*i+1)*Math.pow(v,2*j);}
  for(k=0;k<YT.length;k++){i=YT[k][0];j=YT[k][1];y+=b[k]*Math.pow(u,2*i)*Math.pow(v,2*j+1);}return [x,y];};}
function cylFn(p){var A=p[0],b1=p[1],b2=p[2],b3=p[3];return function(la,lo){var t=la/HP;return [A*lo,b1*t+b2*t*t*t+b3*Math.pow(t,5)];};}
function pseudoFn(p){var a0=p[0],a1=p[1],a2=p[2],b1=p[3],b2=p[4],b3=p[5];return function(la,lo){var t=la/HP,h=a0+a1*t*t+a2*Math.pow(t,4);return [lo*h,b1*t+b2*t*t*t+b3*Math.pow(t,5)];};}
function azimFn(p){var d1=p[0],d2=p[1],d3=p[2];return function(la,lo){var c=HP-la,R=d1*c+d2*c*c+d3*c*c*c;return [R*Math.sin(lo),-R*Math.cos(lo)];};}
function superFn(p){var pw=p[0],A=p[1],B=p[2],b1=p[3],b2=p[4],b3=p[5],nm=(b1+b2+b3)||1e-9;return function(la,lo){var t=la/HP,eta=(b1*t+b2*t*t*t+b3*Math.pow(t,5))/nm;eta=Math.max(-0.999999,Math.min(0.999999,eta));var W=Math.max(1e-9,1-Math.pow(Math.abs(eta),pw));return [(lo/PI)*A*Math.pow(W,1/pw),B*eta];};}
function fromParams(type,params){if(type==='lent')return familyFn(params.slice(0,NX),params.slice(NX));if(type==='cyl')return cylFn(params);if(type==='pseudo')return pseudoFn(params);if(type==='azim')return azimFn(params);if(type==='super')return superFn(params);return CL.Equirectangular;}
function sinc(al){return Math.abs(al)<1e-8?1:al/Math.sin(al);}
function clamp(x,a,b){return x<a?a:(x>b?b:x);}
var CL={
  Equirectangular:function(la,lo){return [lo,la];},
  Mercator:function(la,lo){la=clamp(la,-1.4835,1.4835);return [lo,Math.asinh(Math.tan(la))];},
  Sinusoidal:function(la,lo){return [lo*Math.cos(la),la];},
  Mollweide:function(la,lo){var t=la;if(Math.abs(Math.abs(la)-HP)<1e-9)t=Math.sign(la)*HP;else{for(var n=0;n<20;n++){var f=2*t+Math.sin(2*t)-PI*Math.sin(la),fp=2+2*Math.cos(2*t);if(Math.abs(fp)<1e-12)break;var s=f/fp;t-=s;if(Math.abs(s)<1e-12)break;}}return [(2*SQRT2/PI)*lo*Math.cos(t),SQRT2*Math.sin(t)];},
  Hammer:function(la,lo){var d=Math.sqrt(1+Math.cos(la)*Math.cos(lo/2));return [2*SQRT2*Math.cos(la)*Math.sin(lo/2)/d,SQRT2*Math.sin(la)/d];},
  "Winkel Tripel":function(la,lo){var al=Math.acos(clamp(Math.cos(la)*Math.cos(lo/2),-1,1)),f=sinc(al),p1=Math.acos(2/PI);return [0.5*(lo*Math.cos(p1)+2*Math.cos(la)*Math.sin(lo/2)*f),0.5*(la+Math.sin(la)*f)];},
  "Equal Earth":function(la,lo){var A1=1.340264,A2=-0.081106,A3=0.000893,A4=0.003796,th=Math.asin(clamp(SQRT3/2*Math.sin(la),-1,1)),den=9*A4*Math.pow(th,8)+7*A3*Math.pow(th,6)+3*A2*th*th+A1;return [2*SQRT3*lo*Math.cos(th)/(3*den),A1*th+A2*Math.pow(th,3)+A3*Math.pow(th,7)+A4*Math.pow(th,9)];},
  "Azimuthal Equidist.":function(la,lo){var r=HP-la;return [r*Math.sin(lo),-r*Math.cos(lo)];},
  Robinson:(function(){var X=[1,.9986,.9954,.99,.9822,.973,.96,.9427,.9216,.8962,.8679,.835,.7986,.7597,.7186,.6732,.6213,.5722,.5322],Y=[0,.062,.124,.186,.248,.31,.372,.434,.4958,.5571,.6176,.6769,.7346,.7903,.8435,.8936,.9394,.9761,1];return function(la,lo){var ad=Math.abs(la*180/PI)/5,i=Math.floor(ad),px,py;if(i>=18){px=X[18];py=Y[18];}else{var fr=ad-i;px=X[i]+fr*(X[i+1]-X[i]);py=Y[i]+fr*(Y[i+1]-Y[i]);}return [0.8487*px*lo,1.3523*py*(la>=0?1:-1)];};})()
};
function contenderFn(c){return c.type==='lent'?familyFn(c.a,c.b):(CL[c.classic]||CL.Equirectangular);}

/* ---------- geometry + rendering ---------- */
function buildGeometry(){var G={ocean:[],grat:[],land:[]},b=[],s;
  for(s=-90;s<=90;s+=2)b.push([s,-180]);for(s=-180;s<=180;s+=2)b.push([90,s]);
  for(s=90;s>=-90;s-=2)b.push([s,180]);for(s=-180;s<=180;s+=2)b.push([-90,-s]);G.ocean.push(b);
  for(var lo=-180;lo<=180;lo+=30){var m=[];for(var la=-90;la<=90;la+=3)m.push([la,lo]);G.grat.push(m);}
  for(var la2=-60;la2<=60;la2+=30){var p=[];for(var lo2=-180;lo2<=180;lo2+=3)p.push([la2,lo2]);G.grat.push(p);}
  D.land.forEach(function(r){G.land.push(r.map(function(pt){return [pt[1],pt[0]];}));});return G;}
var GEO=buildGeometry(),GRAT0=1,LAND0=1+GEO.grat.length,ALL=GEO.ocean.concat(GEO.grat,GEO.land);
function fit(fn,w,h,m){var mnx=1e9,mxx=-1e9,mny=1e9,mxy=-1e9;
  for(var la=-90;la<=90;la+=3)for(var lo=-180;lo<=180;lo+=6){var xy=fn(la*PI/180,lo*PI/180);if(xy[0]<mnx)mnx=xy[0];if(xy[0]>mxx)mxx=xy[0];if(xy[1]<mny)mny=xy[1];if(xy[1]>mxy)mxy=xy[1];}
  var sc=Math.min((w-2*m)/(mxx-mnx),(h-2*m)/(mxy-mny));return {sc:sc,ox:(w-sc*(mxx+mnx))/2,oy:(h+sc*(mxy+mny))/2};}
function p2s(fn,ft,la,lo){var xy=fn(la*PI/180,lo*PI/180);return [ft.sc*xy[0]+ft.ox,-ft.sc*xy[1]+ft.oy];}
function screenSet(fn,w,h,m){var ft=fit(fn,w,h,m);return ALL.map(function(poly){return poly.map(function(pt){return p2s(fn,ft,pt[0],pt[1]);});});}
function css(v){return getComputedStyle(document.body).getPropertyValue(v).trim();}
function setupCanvas(cv){var dpr=Math.min(2,window.devicePixelRatio||1),r=cv.getBoundingClientRect();cv.width=Math.max(1,Math.round(r.width*dpr));cv.height=Math.max(1,Math.round(r.height*dpr));var ctx=cv.getContext('2d');ctx.setTransform(dpr,0,0,dpr,0,0);return {ctx:ctx,w:r.width,h:r.height};}
function jump(p,q,w){return Math.abs(p[0]-q[0])>w*0.5;}
function drawSet(ctx,sets,w,h,o){var Co=css('--ocean'),Cl=css('--land'),Cll=css('--landline'),Cli=css('--line');
  ctx.clearRect(0,0,w,h);var ob=sets[0];ctx.beginPath();ob.forEach(function(pt,i){i?ctx.lineTo(pt[0],pt[1]):ctx.moveTo(pt[0],pt[1]);});ctx.closePath();ctx.fillStyle=Co;ctx.fill();
  if(o.grat){ctx.strokeStyle=Cli;ctx.lineWidth=o.mini?0.4:0.7;ctx.globalAlpha=.8;for(var g=GRAT0;g<GRAT0+GEO.grat.length;g++){var pl=sets[g];ctx.beginPath();for(var i=0;i<pl.length;i++){if(i&&jump(pl[i-1],pl[i],w))ctx.moveTo(pl[i][0],pl[i][1]);else i?ctx.lineTo(pl[i][0],pl[i][1]):ctx.moveTo(pl[i][0],pl[i][1]);}ctx.stroke();}ctx.globalAlpha=1;}
  if(o.land!==false){ctx.fillStyle=Cl;ctx.strokeStyle=Cll;ctx.lineWidth=o.mini?0.3:0.5;for(var L=LAND0;L<sets.length;L++){var rg=sets[L];if(rg.length<3)continue;ctx.beginPath();var br=false;for(var j=0;j<rg.length;j++){if(j&&jump(rg[j-1],rg[j],w)){br=true;break;}j?ctx.lineTo(rg[j][0],rg[j][1]):ctx.moveTo(rg[j][0],rg[j][1]);}if(!br){ctx.closePath();ctx.fill();ctx.stroke();}}}}
function mix(a,b,t){var A=hx(a),B=hx(b);return 'rgb('+Math.round(A[0]+(B[0]-A[0])*t)+','+Math.round(A[1]+(B[1]-A[1])*t)+','+Math.round(A[2]+(B[2]-A[2])*t)+')';}
function hx(c){c=c.trim();if(c[0]==='#'){if(c.length===4)c='#'+c[1]+c[1]+c[2]+c[2]+c[3]+c[3];return [parseInt(c.substr(1,2),16),parseInt(c.substr(3,2),16),parseInt(c.substr(5,2),16)];}var m=c.match(/[\d.]+/g);return m?[+m[0],+m[1],+m[2]]:[128,128,128];}
function drawTissot(ctx,fn,ft){var Cs=css('--shape'),Ca=css('--area'),pts=[],mags=[],eh=1e-4,lats=[-60,-30,0,30,60],lons=[-150,-90,-30,30,90,150,180];
  lats.forEach(function(la){lons.forEach(function(lo){var laR=la*PI/180,loR=lo*PI/180,co=Math.cos(laR);if(co<1e-3)return;
    var f=fn(laR,loR),fa=fn(laR+eh,loR),fb=fn(laR-eh,loR),fc=fn(laR,loR+eh),fd=fn(laR,loR-eh);
    var xla=(fa[0]-fb[0])/(2*eh),yla=(fa[1]-fb[1])/(2*eh),xlo=(fc[0]-fd[0])/(2*eh)/co,ylo=(fc[1]-fd[1])/(2*eh)/co;
    var A=[[ft.sc*xla,ft.sc*xlo],[-ft.sc*yla,-ft.sc*ylo]],cx=ft.sc*f[0]+ft.ox,cy=-ft.sc*f[1]+ft.oy;
    var E=A[0][0]*A[0][0]+A[1][0]*A[1][0],Gg=A[0][1]*A[0][1]+A[1][1]*A[1][1],Ff=A[0][0]*A[0][1]+A[1][0]*A[1][1];
    var disc=Math.sqrt(Math.max(0,(E-Gg)*(E-Gg)+4*Ff*Ff)),s1=Math.sqrt(Math.max(0,(E+Gg+disc)/2)),s2=Math.sqrt(Math.max(1e-9,(E+Gg-disc)/2));
    mags.push(Math.sqrt(s1*s2));pts.push({A:A,cx:cx,cy:cy,s1:s1,s2:s2});});});
  if(!pts.length)return;mags.sort(function(a,b){return a-b;});var med=mags[Math.floor(mags.length/2)]||1,rs=Math.min(ft.sc*.2,16)/med;
  pts.forEach(function(p){var om=2*Math.asin(clamp((p.s1-p.s2)/(p.s1+p.s2),0,1)),tt=clamp(om/(45*PI/180),0,1);ctx.beginPath();
    for(var k=0;k<=28;k++){var th=k/28*2*PI,dx=p.A[0][0]*Math.cos(th)+p.A[0][1]*Math.sin(th),dy=p.A[1][0]*Math.cos(th)+p.A[1][1]*Math.sin(th);var X=p.cx+rs*dx,Y=p.cy+rs*dy;k?ctx.lineTo(X,Y):ctx.moveTo(X,Y);}
    ctx.closePath();ctx.fillStyle=mix(Cs,Ca,tt);ctx.globalAlpha=.42;ctx.fill();ctx.globalAlpha=.95;ctx.lineWidth=1;ctx.strokeStyle=mix(Cs,Ca,tt);ctx.stroke();ctx.globalAlpha=1;});}
function drawMini(cv,fn){var s=setupCanvas(cv);drawSet(s.ctx,screenSet(fn,s.w,s.h,7),s.w,s.h,{grat:true,mini:true});}

/* ---------- EXPLORE ---------- */
var famByType={};D.families.forEach(function(f){famByType[f.type]=f;});
var state={type:'lent',w:[1/3,1/3,1/3]};
function idw(fam,w){var pts=fam.points,np=pts[0].params.length,params=new Array(np).fill(0),sh=0,ar=0,di=0,ws=0;
  for(var s=0;s<pts.length;s++){var d2=0;for(var t=0;t<3;t++){var dd=w[t]-pts[s].w[t];d2+=dd*dd;}var wt=1/(Math.pow(d2,3)+1e-9);
    for(var k=0;k<np;k++)params[k]+=wt*pts[s].params[k];sh+=wt*pts[s].shape;ar+=wt*pts[s].area;di+=wt*pts[s].dist;ws+=wt;}
  for(k=0;k<np;k++)params[k]/=ws;return {params:params,shape:sh/ws,area:ar/ws,dist:di/ws};}
var mainCv,mainDim,curSet,tgtSet,curFn,curFt,anim=0,tissot=true;
function setTarget(fn){if(!mainDim)return;tgtSet=screenSet(fn,mainDim.w,mainDim.h,16);curFn=fn;curFt=fit(fn,mainDim.w,mainDim.h,16);if(!curSet)curSet=tgtSet.map(function(pl){return pl.map(function(p){return p.slice();});});if(!anim){anim=1;requestAnimationFrame(step);}}
function step(){if(!curSet||!tgtSet){anim=0;return;}var done=true;for(var i=0;i<curSet.length;i++)for(var j=0;j<curSet[i].length;j++){var dx=tgtSet[i][j][0]-curSet[i][j][0],dy=tgtSet[i][j][1]-curSet[i][j][1];if(Math.abs(dx)>.25||Math.abs(dy)>.25)done=false;curSet[i][j][0]+=dx*.22;curSet[i][j][1]+=dy*.22;}
  renderMain(!done);if(!done)requestAnimationFrame(step);else{anim=0;curSet=tgtSet.map(function(pl){return pl.map(function(p){return p.slice();});});renderMain(false);}}
function renderMain(animating){if(!mainDim)return;drawSet(mainDim.ctx,curSet,mainDim.w,mainDim.h,{grat:true});if(!animating&&tissot&&curFn&&curFt)drawTissot(mainDim.ctx,curFn,curFt);}
function refreshExplore(instant){var fam=famByType[state.type],it=idw(fam,state.w),fn=fromParams(state.type,it.params);
  document.getElementById('mapname').innerHTML=fam.outline.replace(/ \(.*/,'')+' <span class="tag">discovered</span>';
  document.getElementById('mapsub').textContent='shape '+Math.round(state.w[0]*100)+'% · area '+Math.round(state.w[1]*100)+'% · distance '+Math.round(state.w[2]*100)+'%';
  setMetrics(it.shape,it.area,it.dist);
  if(instant&&mainDim){curFn=fn;curFt=fit(fn,mainDim.w,mainDim.h,16);curSet=screenSet(fn,mainDim.w,mainDim.h,16);tgtSet=curSet.map(function(pl){return pl.map(function(p){return p.slice();});});renderMain(false);}else setTarget(fn);
  drawTri();}
function setMetrics(sh,ar,di){var W=D.winkel;
  document.getElementById('mShape').textContent=sh.toFixed(3);document.getElementById('mArea').textContent=ar.toFixed(3);document.getElementById('mDist').textContent=di.toFixed(3);
  document.getElementById('mShapeCmp').innerHTML=cmp((sh-W.shape)/W.shape*100);document.getElementById('mAreaCmp').innerHTML=cmp((ar-W.area)/W.area*100);document.getElementById('mDistCmp').innerHTML=cmp((di-W.dist)/W.dist*100);}
function cmp(p){var c=p<-0.5?'var(--good)':(p>0.5?'var(--bad)':'var(--muted)');return '<b style="color:'+c+'">'+(p<0?'−':'+')+Math.abs(p).toFixed(0)+'%</b> vs WT';}
function setWeights(w,fromTri){var s=w[0]+w[1]+w[2]||1;state.w=[w[0]/s,w[1]/s,w[2]/s];syncSliders();refreshExplore(false);}
function syncSliders(){document.getElementById('slS').value=Math.round(state.w[0]*100);document.getElementById('slA').value=Math.round(state.w[1]*100);document.getElementById('slD').value=Math.round(state.w[2]*100);
  document.getElementById('vS').textContent=Math.round(state.w[0]*100)+'%';document.getElementById('vA').textContent=Math.round(state.w[1]*100)+'%';document.getElementById('vD').textContent=Math.round(state.w[2]*100)+'%';}
/* constrained sliders: moving one rebalances the other two proportionally */
function sliderChange(idx,val){var v=clamp(val/100,0,1),others=[0,1,2].filter(function(i){return i!==idx;});
  var rem=1-v,o0=state.w[others[0]],o1=state.w[others[1]],sum=o0+o1;var nw=[0,0,0];nw[idx]=v;
  if(sum<1e-6){nw[others[0]]=rem/2;nw[others[1]]=rem/2;}else{nw[others[0]]=rem*o0/sum;nw[others[1]]=rem*o1/sum;}
  state.w=nw;syncSliders();refreshExplore(false);}
document.getElementById('slS').addEventListener('input',function(e){sliderChange(0,+e.target.value);});
document.getElementById('slA').addEventListener('input',function(e){sliderChange(1,+e.target.value);});
document.getElementById('slD').addEventListener('input',function(e){sliderChange(2,+e.target.value);});
/* outline chips */
(function(){var box=document.getElementById('famchips');D.families.forEach(function(f){var b=document.createElement('button');b.textContent=f.outline.replace(/ \(.*/,'');b.setAttribute('aria-pressed',f.type===state.type);b.dataset.type=f.type;b.title=f.outline;box.appendChild(b);
  b.addEventListener('click',function(){state.type=f.type;[].forEach.call(box.children,function(c){c.setAttribute('aria-pressed',c.dataset.type===state.type);});refreshExplore(false);});});})();
/* triangle */
var TRI=document.getElementById('tri'),Sc={x:150,y:20},Ac={x:32,y:236},Dc={x:268,y:236};
function bary2xy(w){return {x:w[0]*Sc.x+w[1]*Ac.x+w[2]*Dc.x,y:w[0]*Sc.y+w[1]*Ac.y+w[2]*Dc.y};}
function xy2bary(px,py){var det=(Ac.y-Dc.y)*(Sc.x-Dc.x)+(Dc.x-Ac.x)*(Sc.y-Dc.y);
  var w0=((Ac.y-Dc.y)*(px-Dc.x)+(Dc.x-Ac.x)*(py-Dc.y))/det,w1=((Dc.y-Sc.y)*(px-Dc.x)+(Sc.x-Dc.x)*(py-Dc.y))/det,w2=1-w0-w1;
  return [Math.max(0,w0),Math.max(0,w1),Math.max(0,w2)];}
function drawTri(){var fam=famByType[state.type],Cg=css('--good'),Cf=css('--faint'),Cl=css('--line'),Cs=css('--shape'),Ca=css('--area'),Cd=css('--dist');
  var norms=fam.points.map(function(p){return p.norm;}),nmin=Math.min.apply(null,norms),nmax=Math.max.apply(null,norms),N=14,s='';
  function cell(w){var it=0,ws=0;for(var q=0;q<fam.points.length;q++){var d2=0;for(var t=0;t<3;t++){var dd=w[t]-fam.points[q].w[t];d2+=dd*dd;}var wt=1/(Math.pow(d2,3)+1e-9);it+=wt*fam.points[q].norm;ws+=wt;}return mix(Cg,Cf,clamp((it/ws-nmin)/(nmax-nmin||1),0,1));}
  function pt(i,j){var c=bary2xy([(N-i-j)/N,i/N,j/N]);return c.x+','+c.y;}
  for(var i=0;i<N;i++)for(var j=0;j<N-i;j++){s+='<polygon points="'+pt(i,j)+' '+pt(i+1,j)+' '+pt(i,j+1)+'" fill="'+cell([(N-i-j-.667)/N,(i+.333)/N,(j+.333)/N])+'"/>';if(i+j<N-1)s+='<polygon points="'+pt(i+1,j)+' '+pt(i,j+1)+' '+pt(i+1,j+1)+'" fill="'+cell([(N-i-j-1.333)/N,(i+.667)/N,(j+.667)/N])+'"/>';}
  s+='<polygon points="'+Sc.x+','+Sc.y+' '+Ac.x+','+Ac.y+' '+Dc.x+','+Dc.y+'" fill="none" stroke="'+Cl+'" stroke-width="1.5"/>';
  s+='<text x="'+Sc.x+'" y="'+(Sc.y-6)+'" fill="'+Cs+'" font-size="12" font-weight="600" text-anchor="middle">SHAPE</text>';
  s+='<text x="'+(Ac.x-2)+'" y="'+(Ac.y+16)+'" fill="'+Ca+'" font-size="12" font-weight="600">AREA</text>';
  s+='<text x="'+(Dc.x+2)+'" y="'+(Dc.y+16)+'" fill="'+Cd+'" font-size="12" font-weight="600" text-anchor="end">DISTANCE</text>';
  var h=bary2xy(state.w);s+='<circle cx="'+h.x+'" cy="'+h.y+'" r="7" fill="'+css('--ours')+'" stroke="#fff" stroke-width="2"/>';
  TRI.innerHTML=s;document.getElementById('trihint').textContent='greener = lower total distortion for this outline';}
function triPoint(e){var r=TRI.getBoundingClientRect(),cx=(e.touches?e.touches[0].clientX:e.clientX),cy=(e.touches?e.touches[0].clientY:e.clientY);setWeights(xy2bary((cx-r.left)/r.width*300,(cy-r.top)/r.height*264),true);}
var dragging=false;
TRI.addEventListener('pointerdown',function(e){dragging=true;TRI.setPointerCapture(e.pointerId);triPoint(e);});
TRI.addEventListener('pointermove',function(e){if(dragging)triPoint(e);});
TRI.addEventListener('pointerup',function(){dragging=false;});TRI.addEventListener('pointercancel',function(){dragging=false;});
document.getElementById('tgTissot').addEventListener('change',function(e){tissot=e.target.checked;renderMain(false);});
function resizeExplore(){mainCv=document.getElementById('mainmap');mainDim=setupCanvas(mainCv);curSet=null;refreshExplore(true);}

/* ---------- API client (with localStorage fallback) ---------- */
var api={mode:'?'};
function postVote(w,l){
  return fetch('/api/vote',{method:'POST',headers:{'content-type':'application/json'},body:JSON.stringify({winner:w,loser:l})})
    .then(function(r){if(!r.ok)throw 0;return r.json();}).then(function(){api.mode='live';})
    .catch(function(){api.mode='local';var v=lsVotes();v.push(w+'|'+l);try{localStorage.setItem('pa_votes',JSON.stringify(v.slice(-2000)));}catch(e){}});
}
function lsVotes(){try{return JSON.parse(localStorage.getItem('pa_votes')||'[]');}catch(e){return [];}}
function getLeaderboard(){
  return fetch('/api/leaderboard').then(function(r){if(!r.ok)throw 0;return r.json();}).then(function(j){api.mode='live';return j;})
    .catch(function(){api.mode='local';return {standings:eloFromVotes(lsVotes()),total:lsVotes().length,mode:'local'};});
}
function eloFromVotes(votes){var R={},W={},L={};(D.contenders||[]).forEach(function(c){R[c.key]=1000;W[c.key]=0;L[c.key]=0;});
  votes.forEach(function(s){var p=s.split('|'),w=p[0],l=p[1];if(!(w in R)||!(l in R))return;var ew=1/(1+Math.pow(10,(R[l]-R[w])/400)),K=24;R[w]+=K*(1-ew);R[l]-=K*(1-ew);W[w]++;L[l]++;});
  return (D.contenders||[]).map(function(c){return {key:c.key,elo:R[c.key],wins:W[c.key],losses:L[c.key],games:W[c.key]+L[c.key]};});}
function setStatus(el){el.textContent=api.mode==='live'?'● live · shared':'◦ local · this browser only';el.className='status '+(api.mode==='live'?'live':'local');}

/* ---------- VOTE ---------- */
var CT=D.contenders,byKey={};CT.forEach(function(c){byKey[c.key]=c;});
var pair=[null,null],nVotes=0;
function newPair(){var idx=CT.map(function(_,i){return i;});for(var i=idx.length-1;i>0;i--){var j=Math.floor(Math.random()*(i+1)),t=idx[i];idx[i]=idx[j];idx[j]=t;}
  var a=CT[idx[0]],b=CT[idx[1]];pair=[a,b];
  document.getElementById('nmA').innerHTML=a.name+(a.ours?' <span class="tag2">AI</span>':'');document.getElementById('nmB').innerHTML=b.name+(b.ours?' <span class="tag2">AI</span>':'');
  document.getElementById('whyA').textContent=a.why;document.getElementById('whyB').textContent=b.why;
  drawMini(document.getElementById('canA'),contenderFn(a));drawMini(document.getElementById('canB'),contenderFn(b));}
function vote(i){if(!pair[0])return;var w=pair[i],l=pair[1-i];nVotes++;document.getElementById('voteCount').textContent='you’ve cast '+nVotes+' vote'+(nVotes===1?'':'s')+' this visit';
  postVote(w.key,l.key).then(function(){setStatus(document.getElementById('voteStatus'));});newPair();}
document.getElementById('duelA').addEventListener('click',function(){vote(0);});
document.getElementById('duelB').addEventListener('click',function(){vote(1);});
document.getElementById('duelA').addEventListener('keydown',function(e){if(e.key==='Enter'||e.key===' '){e.preventDefault();vote(0);}});
document.getElementById('duelB').addEventListener('keydown',function(e){if(e.key==='Enter'||e.key===' '){e.preventDefault();vote(1);}});
document.getElementById('skipBtn').addEventListener('click',newPair);
var voteInit=false;
function openVote(){if(!voteInit){voteInit=true;newPair();getLeaderboard().then(function(){setStatus(document.getElementById('voteStatus'));});}else{drawMini(document.getElementById('canA'),contenderFn(pair[0]));drawMini(document.getElementById('canB'),contenderFn(pair[1]));}}

/* ---------- LEADERBOARD ---------- */
function loadBoard(){var grid=document.getElementById('lbgrid');
  getLeaderboard().then(function(j){
    document.getElementById('lbTotal').textContent=(j.total||0)+' vote'+((j.total||0)===1?'':'s')+' so far'+(api.mode==='local'?' (local only — deploy with a store to share).':'.');
    var st=(j.standings||[]).slice().filter(function(r){return byKey[r.key];}).sort(function(a,b){return b.elo-a.elo;});
    var maxE=Math.max.apply(null,st.map(function(r){return r.elo;})),minE=Math.min.apply(null,st.map(function(r){return r.elo;}));
    grid.innerHTML='';
    st.forEach(function(r,i){var c=byKey[r.key];var row=document.createElement('div');row.className='lbrow'+(i===0?' top':'');
      row.innerHTML='<div class="rk">'+(i+1)+'</div><div class="mw"><canvas></canvas></div>'+
        '<div><div class="nm">'+c.name+(c.ours?' <span class="tag2">AI</span>':'')+'</div><div class="meta">'+r.wins+'–'+r.losses+' · '+r.games+' games</div></div>'+
        '<div class="elo">'+Math.round(r.elo)+'<small>Elo</small></div>';
      grid.appendChild(row);drawMini(row.querySelector('canvas'),contenderFn(c));});
  });}

/* ---------- tabs ---------- */
var TABS=['explore','vote','leaderboard'];
function show(tab){TABS.forEach(function(t){document.getElementById(t).hidden=(t!==tab);});
  [].forEach.call(document.querySelectorAll('.tabs button'),function(b){b.setAttribute('aria-selected',b.dataset.tab===tab);});
  if(tab==='explore')resizeExplore();if(tab==='vote')openVote();if(tab==='leaderboard')loadBoard();
  history.replaceState(null,'','#'+tab);window.scrollTo(0,0);}
[].forEach.call(document.querySelectorAll('.tabs button'),function(b){b.addEventListener('click',function(){show(b.dataset.tab);});});
document.getElementById('refreshBtn').addEventListener('click',loadBoard);

/* provenance */
document.getElementById('prov').innerHTML='<b>Pareto Atlas.</b> Maps discovered by pure-Python gradient descent over symmetric-polynomial projection families, scored via Tissot’s indicatrix (shape, area) and great-circle vs. map distance over thousands of point pairs. Coastlines: <code>Natural Earth 110m</code>. Votes are stored server-side (Upstash Redis) when deployed, else in your browser. Source & method: the project README.';

/* boot */
var rt;window.addEventListener('resize',function(){clearTimeout(rt);rt=setTimeout(function(){var cur=TABS.filter(function(t){return !document.getElementById(t).hidden;})[0];if(cur==='explore')resizeExplore();if(cur==='vote'&&pair[0]){drawMini(document.getElementById('canA'),contenderFn(pair[0]));drawMini(document.getElementById('canB'),contenderFn(pair[1]));}if(cur==='leaderboard')loadBoard();},160);});
var start=(location.hash||'#explore').slice(1);show(TABS.indexOf(start)>=0?start:'explore');
})();

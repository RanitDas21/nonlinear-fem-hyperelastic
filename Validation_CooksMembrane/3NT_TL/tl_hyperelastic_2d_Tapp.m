% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function tl_hyperelastic_2d_Tapp()
clear all
%
%
mu = 0.8;
thickness = 1;
Gp = 8000;
%Gp = 5;
%
% -------------------------------------------------------------------------
presnodes = load('fixed_nodes.dat'); 
dispbc = length(presnodes);
dofind = [ones(1,dispbc); ones(1,dispbc)]; 
dispvals = zeros(2,dispbc);
% -------------------------------------------------------------------------
%
node_data = load('nodes.dat');
nnode = size(node_data,1);
%
%
el_data = load('elements.dat');
nelem = size(el_data,1);
%
%
npel = 3; % Nodes per element
ndofn = 2; % No of nodal d-o-fs
ntens = 4; % No of tensor components
ndim = 2; % No of physical dimensions 
ndofel = npel*ndofn; % No of d-o-fs per element
ndoft = ndofn*nnode;
%
% Nodal coordinate data ---------------------------------------------------
% Will be imported from mesh data for 2D and 3D cases ---------------------
ncoord = node_data(:,2:end);
%
% Element connectivity data -----------------------------------------------
% Will be imported from mesh data for 2D and 3D cases ---------------------
elconn = el_data(:,2:end);
%
%
%u1vals = zeros(nnode,1);
%u2vals = zeros(nnode,1);

uinit = zeros(ndoft,1);

% Specify number of Gauss points
ngp=3;
%
%
%
params = [mu thickness Gp];
iparams = [nnode nelem npel ndofn ngp ntens ndim];
%
%
uval=uinit;
Tappval = 0.00:1/320:1/16;
%
%
uext = zeros(length(Tappval),1);
sigma11val = zeros(length(Tappval),1);
%
%
%

for tt=1:1:length(Tappval)
%    
deltau = 100.*ones(ndoft,1); 
Rg = 100.*ones(ndoft,1);
%
incrcount = 0 ; % initialize the count of the number of NR increments to 0.
%
while(max(abs(deltau))>1.0e-6 || max(abs(Rg))>1.0e-8)
%
if incrcount>20
    disp('Too many iterations - NR did not converge.')
    break
end
% Get assembled element force and stiffness matrices ----------------------
[Tg,Rg] = femats(ncoord,elconn,params,iparams,uval);
%
%
% Apply traction boundary conditions  -------------------------------------
% Simply let everything zero if no external point/traction load is aapplied
%
% Consider the total force arising from the traction vector to be simply 
% distributed amongst the boundary nodes associated with the traction
% boundary conditions.
%
trval = Tappval(tt); area = 16;
loadnodes = load('force_nodes.dat'); 
nsnodes = length(loadnodes);
Rg(2*loadnodes(1)) = Rg(2*loadnodes(1)) - 0.5*trval/(nsnodes-1)*area;
Rg(2*loadnodes(end)) = Rg(2*loadnodes(end)) - 0.5*trval/(nsnodes-1)*area;
%
for jk = 2:1:nsnodes-1
  Rg(2*loadnodes(jk)) = Rg(2*loadnodes(jk)) - trval/(nsnodes-1)*area;
end
%
% -------------------------------------------------------------------------
% Apply displacement boundary conditions ----------------------------------
% -------------------------------------------------------------------------
%
for ii = 1:1:dispbc
    for jj=1:1:ndofn  
    if dofind(jj,ii)==1    
       Tg(ndofn*(presnodes(ii)-1)+jj,:) = zeros(1,ndofn*nnode);
       Rg(ndofn*(presnodes(ii)-1)+jj) = uval(ndofn*(presnodes(ii)-1)+jj) - dispvals(jj,ii);
    end
    end
end
%
for ii = 1:1:dispbc
    for jj=1:1:ndofn  
    if dofind(jj,ii)==1    
       Tg(ndofn*(presnodes(ii)-1)+jj,ndofn*(presnodes(ii)-1)+jj) = 1;
     end
    end
end
%
%
deltau = - Tg\Rg;
%
%
uval = uval+deltau;
incrcount = incrcount +1;
%
end

if incrcount>20
    disp('Programme aborted.')
    break
end

disp(strcat('increment number: ',  num2str(tt)));


plotdata = plotop(ncoord,elconn,params,iparams,uval,'N','meshdef') ;

dispu2(tt) = uval(2*3);


end
%
% Postprocessing ----------------------------------------------------------
%
%
% 
plotdata = plotop(ncoord,elconn,params,iparams,uval,'Y','sigma12') ;

%

figure

plot(Tappval,dispu2,'m-')
%
%
%
%
%
% -------------------------------------------------------------------------
end
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%
%
% Global force and stiffness matrix calculation subroutine ----------------
%
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%
function [Tg,Rg] = femats(ncoord,elconn,params,iparams,uval)
%
nnode = iparams(1);
nelem = iparams(2);
npel = iparams(3);
ndofn = iparams(4);
%
Tg = zeros(nnode*ndofn,nnode*ndofn);
Rg = zeros(nnode*ndofn,1);
%
for j=1:1:nelem
coords = [];
uvalvec = [];

for ii=1:1:npel
coords = [coords ; ncoord(elconn(j,ii),:)];
uvalvec = [uvalvec; uval(ndofn*elconn(j,ii)-1,1); uval(ndofn*elconn(j,ii),1)];
end
%
%
% Call function for element stiffness and force matrices ------------------
[telem,relem] = kelmats(coords,params,iparams,uvalvec);
%
% Assemble element stiffness and force matrices ---------------------------
%
for ii=1:1:npel
    %
    for idf=1:1:ndofn
        Rg(ndofn*(elconn(j,ii)-1)+idf,1) = Rg(ndofn*(elconn(j,ii)-1)+idf,1) ...
                                           + relem(ndofn*(ii-1)+idf,1);
    end
    %
    for jj=1:1:npel
    for idf=1:1:ndofn
        for jdf =1:1:ndofn
           Tg(ndofn*(elconn(j,ii)-1)+idf,ndofn*(elconn(j,jj)-1)+jdf) = ...
           Tg(ndofn*(elconn(j,ii)-1)+idf,ndofn*(elconn(j,jj)-1)+jdf) + ...
           telem(ndofn*(ii-1)+idf,ndofn*(jj-1)+jdf);    
        end
    end
    end
    %
end
%
end

end
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%
%
% Element force and stiffness matrix calculation subroutine ---------------
%
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%

function [telem,relem] = kelmats(coords,params,iparams,uvalvec)
%
nnode = iparams(1);
nelem = iparams(2);
npel = iparams(3);
ndofn = iparams(4);
ngp = iparams(5);
ntens = iparams(6);
ndim = iparams(7);
%
ndofel = ndofn*npel;
%
% Material parameters -----------------------------------------------------
%
telem = zeros(ndofel,ndofel);
relem = zeros(ndofel,1);
%
% Obtain gauss points and weightages from the kgaussp function ------------
mu = params(1);
thickness = params(2);
Gp = params(3);

[xgp,egp,wgp]=kgaussp(ngp);
%
for jj=1:1:ngp
%
xig = xgp(jj); etg = egp(jj);
%
N = zeros(1,npel); dNdxi = zeros(1,npel); dNdet = zeros(1,npel);
%
% Linear shape functions in natural coordinates ---------------------------
N(1) = 1-xig-etg; N(2) = xig; N(3) = etg;
dNdxi(1) = -1; dNdxi(2) = 1; dNdxi(3) = 0;
dNdet(1) = -1; dNdet(2) = 0; dNdet(3) = 1;
% -------------------------------------------------------------------------
%
% Construct the Jacobian matrix
x1p = 0; dx1pdxi = 0; dx1pdet = 0;
x2p = 0; dx2pdxi = 0; dx2pdet = 0;
%
for ii=1:1:npel
%    
x1p = x1p+N(ii)*coords(ii,1);
dx1pdxi = dx1pdxi + dNdxi(ii)*coords(ii,1);
dx1pdet = dx1pdet + dNdet(ii)*coords(ii,1);
%
x2p = x2p+N(ii)*coords(ii,2);
dx2pdxi = dx2pdxi + dNdxi(ii)*coords(ii,2);
dx2pdet = dx2pdet + dNdet(ii)*coords(ii,2);
%
end
Jac = [dx1pdxi dx2pdxi; dx1pdet dx2pdet];
%
%
% Construct the JI matrix -------------------------------------------------
Jinv = inv(Jac);
%
JI = [Jinv(1,1) Jinv(1,2) 0 0; 
      Jinv(2,1) Jinv(2,2) 0 0; 
      0 0 Jinv(1,1) Jinv(1,2);
      0 0 Jinv(2,1) Jinv(2,2)];
% Construct the NG matrix -------------------------------------------------
NG = [dNdxi(1) 0 dNdxi(2) 0 dNdxi(3) 0;
      dNdet(1) 0 dNdet(2) 0 dNdet(3) 0;
      0 dNdxi(1) 0 dNdxi(2) 0 dNdxi(3);
      0 dNdet(1) 0 dNdet(2) 0 dNdet(3)]; % For linear shape functions
%
% Construct the G matrix --------------------------------------------------
G = JI*NG;
%
% Construct the S and DSDF matrices ---------------------------------------
Gradu_vec = G*uvalvec;
%
Ftens = [1+Gradu_vec(1) Gradu_vec(2) 0; Gradu_vec(3) Gradu_vec(4)+1 0; 0 0 1];
Idtens = eye(3,3);
Jval = det(Ftens);
Finv = inv(Ftens);
%
Stens = mu.*Ftens-mu.*Finv'+Gp*Jval*(Jval-1).*Finv'; % PK-1 for Neo-Hookean
Svc = [Stens(1,1); Stens(1,2); Stens(2,1); Stens(2,2)];
%
rmgp = det(Jac).*G'*Svc;
%
DSDF = zeros(ndim,ndim,ndim,ndim);
%
for iv=1:1:ndim
    for jv=1:1:ndim
        for kv = 1:1:ndim
            for lv = 1:1:ndim
                DSDF(iv,jv,kv,lv) = mu*Idtens(iv,kv)*Idtens(jv,lv) + ...
                                    mu*Finv(lv,iv)*Finv(jv,kv) + ...
                                    Gp*(2*Jval-1)*Finv(lv,kv)*Finv(jv,iv) - ...
                                    Gp*Jval*(Jval-1)*Finv(lv,iv)*Finv(jv,kv);
            end
        end
    end
end

Dmat = [DSDF(1,1,1,1) DSDF(1,1,1,2) DSDF(1,1,2,1) DSDF(1,1,2,2);
           DSDF(1,2,1,1) DSDF(1,2,1,2) DSDF(1,2,2,1) DSDF(1,2,2,2);
           DSDF(2,1,1,1) DSDF(2,1,1,2) DSDF(2,1,2,1) DSDF(2,1,2,2);
           DSDF(2,2,1,1) DSDF(2,2,1,2) DSDF(2,2,2,1) DSDF(2,2,2,2)];

tmgp = det(Jac).*G'*Dmat*G;

% Sumation for performing the Gauss integration ---------------------------
telem = telem + wgp(jj)*tmgp;
relem = relem + wgp(jj)*rmgp;
%
end
%
%
end
%
% -------------------------------------------------------------------------
% Gauss integration points ------------------------------------------------
% -------------------------------------------------------------------------
%
function [xgp,egp,wgp]=kgaussp(ngp)
%
xgp = zeros(1,ngp); wgp = zeros(1,ngp); egp = zeros(1,ngp);
%
if ngp == 1
 xgp = 1/3; egp = 1/3; wgp = 1/2;
elseif ngp==3
 xgp = [1/6 2/3 1/6]; egp=[1/6 1/6 2/3]; wgp = [1/6 1/6 1/6];
else
 disp('Error -- Gauss points not supplied')
 pause
end
%
end
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%
%
% Postprocessing subroutines starts ---------------------------------------
%
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%
%
%
function plotdata = plotop(ncoord,elconn,params,iparams,uout,pltind,sigstr)
%
nnode = iparams(1);
nelem = iparams(2);
npel = iparams(3);
ndofn = iparams(4);
ngp = iparams(5);
%
u1vals = uout(1:2:end);
u2vals = uout(2:2:end);
%
% Calculate displaced coordinates to plot the deformed configuration ------
ncoordd(:,1) = ncoord(:,1)+u1vals;
ncoordd(:,2) = ncoord(:,2)+u2vals;
% -------------------------------------------------------------------------
sigvals = extractstress(ncoord,elconn,params,iparams,uout);
sigvnodes = zeros(nnode,3);


%
%xbound = [];
%sigbound11 = []; sigbound22 = [];




%Find nodal sigma values as average of stress in the elememts it shares --

for jj=1:1:nnode
[posn,~] = find(elconn==jj);
sig11 = 0; sig22 = 0; sig12 = 0;
for ii=1:1:length(posn)
    sig11 = sig11+sigvals(posn(ii),3);
    sig22 = sig22+sigvals(posn(ii),4);
    sig12 = sig12+sigvals(posn(ii),5);
end
sig11 = sig11/length(posn);
sig22 = sig22/length(posn);
sig12 = sig12/length(posn);
sigvnodes(jj,:) = [sig11 sig22 sig12];

% if abs(ncoord(jj,1)-0.0)<1.0e-8 && ncoord(jj,2)>0
% xbound = [xbound; ncoord(jj,2)];
% sigbound11 = [sigbound11; sig11];
% sigbound22 = [sigbound22; sig22];
% end

end

% -------------------------------------------------------------------------
%
%
%
xx00 = [];
xxmod = [];
f = [];
%
map_sigma11 = [];
map_sigma22 = [];
map_sigma12 = [];

%
for jj=1:1:nelem
coords = []; dcoords = [];
for ii=1:1:npel
coords = [coords ; ncoord(elconn(jj,ii),:)];
dcoords = [dcoords ; ncoordd(elconn(jj,ii),:)];
end
%
x10 = coords(:,1);
x20 = coords(:,2);
xx00 = [xx00 ; x10 x20];
%
% 
x10 = dcoords(:,1);
x20 = dcoords(:,2);
%
xxmod = [xxmod ; x10 x20];
fmt = zeros(1,npel);
%
%
for ii=1:1:npel
fmt(1,ii) = npel*jj-(npel-ii);    
map_sigma11 = [map_sigma11; sigvnodes(elconn(jj,ii),1)];
map_sigma22 = [map_sigma22; sigvnodes(elconn(jj,ii),2)];
map_sigma12 = [map_sigma12; sigvnodes(elconn(jj,ii),3)];
end
f = [f; fmt];
%
%
end
%

%
if pltind == 'Y'
% if sigstr == 'sigma11'
%     patch('Faces',f,'Vertices',xxmod,'FaceVertexCData',map_sigma11,...
%     'FaceColor','interp','EdgeColor','none')
%     title('Contours of $\sigma_{11}$ (zoomed)','Interpreter','latex')
%     xlim([-0.2,0.2]);
%     ylim([-0.2,0.2]);
% elseif sigstr == 'sigma22'
%      patch('Faces',f,'Vertices',xxmod,'FaceVertexCData',map_sigma22,...
%     'FaceColor','interp','EdgeColor','none')  
%      title('Contours of $\sigma_{22}$ (zoomed)','Interpreter','latex')
%     xlim([-0.2,0.2]);
%     ylim([-0.2,0.2]);
% elseif sigstr == 'sigma12'
%      patch('Faces',f,'Vertices',xxmod,'FaceVertexCData',map_sigma12,...
%     'FaceColor','interp','EdgeColor','none')  
%      title('Contours of $\sigma_{12}$ (zoomed)','Interpreter','latex')
%     xlim([-0.2,0.2]);
%     ylim([-0.2,0.2]);
% elseif sigstr == 'meshdef'
%     patch('Faces',f,'Vertices',xxmod,...
%        'FaceColor','none','EdgeColor','blue')      
%     hold on
%     patch('Faces',f,'Vertices',xx00,...
%        'FaceColor','none','EdgeColor','black')      
%     hold off
% end
% box on
% axis equal
% colormap jet
%
%
%figure
if sigstr == 'sigma11'
    patch('Faces',f,'Vertices',xxmod,'FaceVertexCData',map_sigma11,...
    'FaceColor','interp','EdgeColor','none')
    title('Contours of $\sigma_{11}$ ','Interpreter','latex')
elseif sigstr == 'sigma22'
     patch('Faces',f,'Vertices',xxmod,'FaceVertexCData',map_sigma22,...
    'FaceColor','interp','EdgeColor','none')   
     title('Contours of $\sigma_{22}$ ','Interpreter','latex')
elseif sigstr == 'sigma12'
     patch('Faces',f,'Vertices',xxmod,'FaceVertexCData',map_sigma12,...
    'FaceColor','interp','EdgeColor','none') 
    title('Contours of $\sigma_{12}$ ','Interpreter','latex') 
elseif sigstr == 'meshdef'
    qq = patch('Faces',f,'Vertices',xxmod,...
       'FaceColor','none','EdgeColor','blue') ;     
    hold on
    patch('Faces',f,'Vertices',xx00,...
       'FaceColor','none','EdgeColor','black')   
    xlim([0,3]);
    hold off
end
%
%

box on
axis equal
colormap jet


end


plotdata=sigvals;


end
%
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%
%
function sigvals = extractstress(ncoord,elconn,params,iparams,uvals)
%
nnode = iparams(1);
nelem = iparams(2);
npel = iparams(3);
ndofn = iparams(4);
ngp = iparams(5);
%
sigvals = [];
for j=1:1:nelem
coords = []; uvel = [];
for ii=1:1:npel
coords = [coords ; ncoord(elconn(j,ii),:)];
uvel = [uvel; uvals(2*elconn(j,ii)-1) ; uvals(2*elconn(j,ii))];
end

selem = sigcalc(coords,params,iparams,uvel);

sigvals = [sigvals; selem];

end

end
% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------
function selem = sigcalc(coords,params,iparams,uvel)
%
%
nnode = iparams(1);
nelem = iparams(2);
npel = iparams(3);
ndofn = iparams(4);
ngp = iparams(5);
ndim = iparams(7);
%
%
ndofel = ndofn*npel;
%
%
% Material parameters -----------------------------------------------------
%
mu = params(1);
thickness = params(2);
Gp = params(3);
%
%
% Obtain gauss points and weightages from the kgaussp function ------------
[xgp,egp,wgp]=kgaussp(ngp);
%
selem = zeros(1,5);
for jj=1:1:ngp
%
xig = xgp(jj); etg = egp(jj);
%
N = zeros(1,npel); dNdxi = zeros(1,npel); dNdet = zeros(1,npel);
%
% Linear shape functions in natural coordinates ---------------------------
N(1) = 1-xig-etg; N(2) = xig; N(3) = etg;
dNdxi(1) = -1; dNdxi(2) = 1; dNdxi(3) = 0;
dNdet(1) = -1; dNdet(2) = 0; dNdet(3) = 1;
% -------------------------------------------------------------------------
%
% Construct the Jacobian matrix
x1p = 0; dx1pdxi = 0; dx1pdet = 0;
x2p = 0; dx2pdxi = 0; dx2pdet = 0;
%
for ii=1:1:npel
%    
x1p = x1p+N(ii)*coords(ii,1);
dx1pdxi = dx1pdxi + dNdxi(ii)*coords(ii,1);
dx1pdet = dx1pdet + dNdet(ii)*coords(ii,1);
%
x2p = x2p+N(ii)*coords(ii,2);
dx2pdxi = dx2pdxi + dNdxi(ii)*coords(ii,2);
dx2pdet = dx2pdet + dNdet(ii)*coords(ii,2);
%
end
Jac = [dx1pdxi dx2pdxi; dx1pdet dx2pdet];
%
%
% Construct the JI matrix -------------------------------------------------
Jinv = inv(Jac);
%
JI = [Jinv(1,1) Jinv(1,2) 0 0; 
      Jinv(2,1) Jinv(2,2) 0 0; 
      0 0 Jinv(1,1) Jinv(1,2);
      0 0 Jinv(2,1) Jinv(2,2)];
% Construct the NG matrix -------------------------------------------------
NG = [dNdxi(1) 0 dNdxi(2) 0 dNdxi(3) 0;
      dNdet(1) 0 dNdet(2) 0 dNdet(3) 0;
      0 dNdxi(1) 0 dNdxi(2) 0 dNdxi(3);
      0 dNdet(1) 0 dNdet(2) 0 dNdet(3)]; % For linear shape functions
%
% Construct the G matrix --------------------------------------------------
G = JI*NG;
% Construct the Kp matrix -------------------------------------------------
Gradu_vec = G*uvel;
%
Ftens = [1+Gradu_vec(1) Gradu_vec(2) 0; Gradu_vec(3) Gradu_vec(4)+1 0; 0 0 1];
Idtens = eye(3,3);
Jval = det(Ftens);
Finv = inv(Ftens);
%
Stens = mu.*Ftens-mu.*Finv'+Gp*Jval*(Jval-1).*Finv';
sigma = (1/Jval).*Stens*Ftens';
%
selem = selem + [x1p x2p sigma(1,1) sigma(2,2) sigma(1,2)];
%
end
%
selem = selem./ngp;
%
%
end
%
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%
% Postprocessing ends -----------------------------------------------------
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%
%
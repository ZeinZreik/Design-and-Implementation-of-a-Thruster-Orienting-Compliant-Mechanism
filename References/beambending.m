function [Fx, Fy, varargout] = beambending(dx, dy, l, h, w, E, varargin)
% This function calculate the bending of a fixed-guided beam using elliptic
% integral solutions. Usage:
%
% [Fx, Fy] = beambending(dx, dy, l, h, w, E)
% This will calculate the x- and y- direction forces which will create
% deflections of dx and dy, where dx and dy may be vectors of x- and y-
% displacements. If dx and dy are vectors, they must be the same size, and
% Fx and Fy will be the same size. The beam is of length l, with a 
% rectangular cross-section and Young's modulus E. The cross-section 
% dimension in the plane of bending is h, and the other cross-section 
% dimension is w.
%
% [Fx, Fy, variable outputs] = beambending(dx, dy, l, h, w, E, npts)
% Additional outputs may also be desired. These outputs will be returned in
% the following order, with the number of outputs returned determined by
% the number of outputs in the calling function:
% errx - the error in approximating the deflection dx. Since the
%   function uses nonlinear equation solutions to find the beam deflections
%   corresponding to dx and dy, errx tells how well each deflection is
%   approximated in the x-direction. Will be the same size as dx.
% erry - the error in approximating the deflection dy. Since the
%   function uses nonlinear equation solutions to find the beam deflections
%   corresponding to dx and dy, erry tells how well each deflection is
%   approximated in the y-direction.
% s - coordinates along the beam length. If dx and dy are scalars, s will
%   be a column vector of length npts. If dx and dy are vectors of length n,
%   s will be a matrix with npts rows and n columns, where each column
%   corresponds to an element of dx and dy.
% M - internal moment in the beam. M will be the same size as s, and gives
%   the internal moment at each beam coordinate. The first and last entries
%   of each column give the reaction moments at each end of the beam.
% Sax - axial stress in the beam. Sax will be the same size as s.
% Sbend - bending stress in the beam. Sbend will be the same size as s.
% x - x-coordinate corresponding to each coordinate along the beam, 
%   assuming the origin is at the fixed end of the beam. x will be the same
%   size as s.
% y - y-coordinate corresponding to each coordinate along the beam. y will
%   be the same size as s. Taken together, x and y give the deflected beam
%   shape.
% mode - deflected beam mode shape. Gives the mode of the beam (1 or 2) for
%   each deflected point. Will be the same size as dx and dy.
% 

tol = 1e-4;
slen = h^2/(12*l^2);
if nargin==6
    npts = 50;
else
    npts = varargin{1};
end
b = abs(dy/l);
a = 1+dx/l;
if a(1)<1
    guess = [180/360 0.05 1];
else
    guess = [5/360 0.9995 1];
end
psi = zeros(size(b));
alpha = psi;
mode = psi;
k1 = psi;
aout = psi;
bout = psi;
aerr = psi;
berr = psi;
sax = zeros(npts,length(b));
sbend = sax;
x = sax;
y = x;
s = y;
beta = y;

for i = 1:length(b)
    mode(i) = findmode(a(i),b(i),slen,npts);
    %if guess(1) > 90/360
        [psi(i),k1(i),alpha(i),aout(i),bout(i),aerr(i),berr(i),sax(:,i),sbend(:,i),x(:,i),y(:,i),s(:,i),beta(:,i)] = solvebeambending(a(i),b(i),mode(i),slen,guess(1:2),npts);
    %else
    toterror = sqrt((a(i)-aout(i))^2 + (b(i)-bout(i))^2);
    if toterror > .0001
        [psi(i),k1(i),alpha(i),aout(i),bout(i),aerr(i),berr(i),sax(:,i),sbend(:,i),x(:,i),y(:,i),s(:,i),beta(:,i)] = solvebeambending2(a(i),b(i),slen,[guess(1) guess(3)],npts);
        %[a(i)*l b(i)*l toterror i]
    end
    guess = [psi(i)/360 k1(i) alpha(i)];
end
sax = E*sax;
sbend = E*sbend;
I = w*h^3/12;
F = alpha*E*I/l^2;
Fx = F.*cosd(psi);
Fy = F.*sind(psi);
s = s*l;
M = beta*E*I/l;
x = x*l;
y = y*l;
errx = (aout - a)*l;
erry = (bout - b)*l;
aerr = aerr*l./sqrt(dx.^2 + dy.^2);
berr = berr*l./sqrt(dx.^2 + dy.^2);
if max(abs(aerr))>0.1 || max(abs(berr))>0.1
    display('Warning: Displacements due to axial stretching or compression are more than 10% of total displacments.')
end
oppos = find(dy<0);
if ~isempty(oppos)
    Fy(oppos) = -Fy(oppos);
    erry(oppos) = -erry(oppos);
    M(oppos) = -M(oppos);
    sbend(:,oppos) = -sbend(:,oppos);
    y(:,oppos) = -y(:,oppos);
    berr(oppos) = -berr(oppos);
end
outvals{1} = errx;
outvals{2} = erry;
outvals{3} = s;
outvals{4} = M;
outvals{5} = sax;
outvals{6} = sbend;
outvals{7} = x;
outvals{8} = y;
outvals{9} = mode;

outvals{10} = k1;
outvals{11} = aerr;
outvals{12} = berr;

varargout = {outvals{1:nargout-2}};

function mode = findmode(ades,bdes,slen,npts)

amax = 1-4*pi^2*slen; % Find the largest value a can take in mode 2
if ades>amax
    mode = 1;
else
    psi = fzero(@findpsi,[180 270],[],ades,slen,npts);
    k = sqrt((1+cosd(psi))/2);

    phi1 = -pi/2;
    phi2 = phi1+2*pi;

    [f1,e1] = elliptic12(phi1,k.^2);
    [f2,e2] = elliptic12(phi2,k.^2);
    alpha = (f2 - f1).^2;
    a = -(cosd(psi).*(f1 - f2 + 2*(e2 - e1)) + ...
        2*sind(psi).*k.*(cos(phi2) - cos(phi1)))./sqrt(alpha);
    b = -(sind(psi).*(f1 - f2 + 2*(e2 - e1)) + ...
        2*k.*cosd(psi).*(cos(phi1)-cos(phi2)))./sqrt(alpha);
    [x,y,s,beta,theta] = guidedbeamshape(psi,phi1,phi2,k,alpha,f1,e1,npts);
    sa = cosd(psi)*(cos(theta)).^2 + sind(psi)*sin(theta).*cos(theta);
    sb = cosd(psi)*cos(theta).*sin(theta) + sind(psi)*(sin(theta)).^2;

    da = slen*alpha.*trapz(s,sa);
    db = slen*alpha.*trapz(s,sb);
    aout = a + da(end);
    bout = b + db(end);
    if bdes>bout
        mode = 1;
    else
        mode = 2;
    end
end

function errout = findpsi(psi,ades,slen,npts)
k = sqrt((1+cosd(psi))/2);

phi1 = -pi/2;
phi2 = phi1+2*pi;

[f1,e1] = elliptic12(phi1,k.^2);
[f2,e2] = elliptic12(phi2,k.^2);
alpha = (f2 - f1).^2;
a = -(cosd(psi).*(f1 - f2 + 2*(e2 - e1)) + ...
    2*sind(psi).*k.*(cos(phi2) - cos(phi1)))./sqrt(alpha);
[x,y,s,beta,theta] = guidedbeamshape(psi,phi1,phi2,k,alpha,f1,e1,npts);
sa = cosd(psi)*(cos(theta)).^2 + sind(psi)*sin(theta).*cos(theta);

da = slen*alpha.*trapz(s,sa);
aout = a + da(end);
errout = aout - ades;

function [psi,k,alpha,aout,bout,aerr,berr,sigmaten,sigmabend,x,y,s,beta] = solvebeambending(a,b,mode,slen,guess,npts)
options = optimset('TolFun',1e-9,'MaxFunEvals',10000,'MaxIter',10000,'TolX',1e-9,'Display','off');
[outvals] = fsolve(@solvefn,guess,options,a,b,mode,slen,npts);
%values
[as,bs,alpha,s,theta,x,y,beta] = guidednonlin(outvals(2),outvals(1)*360,mode,npts);
psi = outvals(1)*360;
k = outvals(2);
if x==0
    aout = as;
    bout = bs;
    aerr = 0;
    berr = 0;
    sigmaten = zeros(size(theta));
    sigmabend = zeros(size(theta));
    beta = zeros(size(theta));
else
    sigmaten = slen*alpha.*cos(psi*pi/180 - theta);
    sigmabend = sqrt(3*slen)*beta;
    sa = cosd(psi)*(cos(theta)).^2 + sind(psi)*sin(theta).*cos(theta);
    sb = cosd(psi)*cos(theta).*sin(theta) + sind(psi)*(sin(theta)).^2;
    da = slen*alpha.*cumtrapz(s,sa);
    db = slen*alpha.*cumtrapz(s,sb);
    aerr = da(end);
    berr = db(end);
    aout = as + da(end);
    bout = bs + db(end);
    x = x + da;
    y = y + db;
end

function [psi,k,alpha,aout,bout,aerr,berr,sigmaten,sigmabend,x,y,s,beta] = solvebeambending2(a,b,slen,guess,npts)
guess(1) = guess(1)*2*pi;
[psi,k,mode] = fixedguidedforceinputsmall(guess,a,b,slen);
psi = psi*180/pi;
[as,bs,alpha,s,theta,x,y,beta] = guidednonlin(k,psi,mode,npts);

if x==0
    aout = as;
    bout = bs;
    aerr = 0;
    berr = 0;
    sigmaten = zeros(size(theta));
    sigmabend = zeros(size(theta));
    beta = zeros(size(theta));
else
    sigmaten = slen*alpha.*cos(psi*pi/180 - theta);
    sigmabend = sqrt(3*slen)*beta;
    sa = cosd(psi)*(cos(theta)).^2 + sind(psi)*sin(theta).*cos(theta);
    sb = cosd(psi)*cos(theta).*sin(theta) + sind(psi)*(sin(theta)).^2;
    da = slen*alpha.*cumtrapz(s,sa);
    db = slen*alpha.*cumtrapz(s,sb);
    aerr = da(end);
    berr = db(end);
    aout = as + da(end);
    bout = bs + db(end);
    x = x + da;
    y = y + db;
end

function about = solvefn(x,ain,bin,mode,slen,npts)
tol = 1e-4;
psi = x(1)*360;
k = x(2);
kmin = sqrt((cosd(psi) + 1)/2);
if k<=kmin
    newa = 1+(kmin-k);
    newb = 0-(kmin-k);
elseif k>=1
    newa = 1+(k-1);
    newb = -(k-1);
else
    [a,b,alpha,s,theta] = guidednonlin(k,psi,mode,npts);
    sa = cosd(psi)*(cos(theta)).^2 + sind(psi)*sin(theta).*cos(theta);
    sb = cosd(psi)*cos(theta).*sin(theta) + sind(psi)*(sin(theta)).^2;
    isa = (sa(1:end-1,:) + sa(2:end,:)).*(s(2:end,:) - s(1:end-1,:))/2;
    isb = (sb(1:end-1,:) + sb(2:end,:)).*(s(2:end,:) - s(1:end-1,:))/2;
    da = slen*alpha.*sum(isa,1);
    db = slen*alpha.*sum(isb,1);
    newa = a + da;
    newb = b + db;
end
if abs(ain-1)>tol
    about(1) = (ain - newa)/(ain-1);
else
    about(1) = (ain - newa)/tol;
end
about(2) = (bin - newb)/bin;


function [a,b,alpha,s,theta,x,y,beta] = guidednonlin(k,psi,mode,npts)
%This function calculates end coordinates and loads for a beam.
%Inputs are k, the elliptic integral modulus, and psi, the angle of the
%applied force, and mode, an integer giving the bending mode. npts is the
%number of points to use in giving the beam shape and moment.
%Outputs are a, the nondimensional horizontal displacement, b, the
%nondimensional vertical displacment, and alpha, the nondimensional load.
%beta1 and beta2 are the nondimensional moment at the fixed and guided beam
%ends, respectively.
%Additionally, outputs x and y give the beam shape.
%Assumptions are an initially-straight guided beam with end loads and moments.

kmin = sqrt((cosd(psi) + 1)/2);
if k<=kmin
    a = 1+10*(kmin-k);
    b = 0-10*(kmin-k);
    alpha = 0;
    x = zeros(npts,1);
    y = x;
    s = x;
    beta = x;
    theta = x;
    return
elseif k>=1
    a = 1+10*(k-1);
    b = -10*(k-1);
    alpha = 0;
    x = zeros(npts,1);
    y = x;
    s = x;
    beta = x;
    theta = x;
    return
end
if sind(psi)>=0
    phi1 = asin(sqrt((cosd(psi) + 1)./(2*k.^2)));
else
    phi1 = -asin(sqrt((cosd(psi) + 1)./(2*k.^2)));
end
switch mode
    case 1
        phi2 = pi - phi1;%first mode
    case 2
        phi2 = 2*pi+phi1;%second mode
    case 3
        phi2 = 3*pi - phi1;%third mode
    case 4
        phi2 = 4*pi+phi1;%fourth mode
    case 5
        phi2 = 5*pi - phi1;%fifth mode
    case 6
        phi2 = 6*pi + phi1;%sixth mode
    case 7
        phi2 = 7*pi - phi1;%seventh mode
    case 8
        phi2 = 8*pi + phi1;%eigth mode
    case 9 %third mode, second type
        phi = asin(sqrt((cosd(psi) + 1)./(2*k.^2)));
        if sind(psi)>=0            
            phi1 = pi-phi;
            phi2 = 4*pi + phi;
        else
            phi1 = pi + phi;
            phi2 = 4*pi - phi;
        end
end
[f1,e1] = elliptic12(phi1,k.^2);
[f2,e2] = elliptic12(phi2,k.^2);
alpha = (f2 - f1).^2;
a = -(cosd(psi).*(f1 - f2 + 2*(e2 - e1)) + ...
    2*sind(psi).*k.*(cos(phi2) - cos(phi1)))./sqrt(alpha);
b = -(sind(psi).*(f1 - f2 + 2*(e2 - e1)) + ...
    2*k.*cosd(psi).*(cos(phi1)-cos(phi2)))./sqrt(alpha);
[x,y,s,beta,theta] = guidedbeamshape(psi,phi1,phi2,k,alpha,f1,e1,npts);


function [x,y,s,beta,theta] = guidedbeamshape(psi,phi1,phi2,k,alpha,f1,e1,npts)

if size(phi1,1)>1
    psi = psi';
    phi1 = phi1';
    phi2 = phi2';
    k = k';
    alpha = alpha';
    f1 = f1';
    e1 = e1';
end

numpts = npts-1;
phi = zeros(numpts+1,length(phi1));
for counter = 1:length(phi1)
phi(:,counter) = (phi1(counter):(phi2(counter) - phi1(counter))/numpts:phi2(counter))';
end
temp2 = ones(size(phi,1),1);
if length(psi)==1
    psi = psi*ones(size(phi1));
end
[f2,e2] = elliptic12(phi,temp2*(k.^2));
x = -(temp2*cosd(psi).*(2*(e2 - temp2*e1)-f2 + temp2*f1) ...
    + 2*(temp2*sind(psi)).*(temp2*k).*(cos(phi) - cos(temp2*phi1)))./sqrt(temp2*alpha);
y = -((temp2*sind(psi)).*(2*(e2 - temp2*e1)-f2 + temp2*f1) - ...
    2*(temp2*k).*(temp2*cosd(psi)).*(cos(phi) - cos(temp2*phi1)))./sqrt(temp2*alpha);
s = (f2 - temp2*f1)./sqrt(temp2*alpha);
beta = 2*(temp2*k).*(temp2*sqrt(alpha)).*cos(phi);
theta = 2*(asin((temp2*k).*sin(phi))-asin((temp2*k).*sin(temp2*phi1)));

function [psi,k,mode] = fixedguidedforceinputsmall(guess,a,b,slen)
options = optimset('MaxFunEvals',10000,'MaxIter',10000,'Display','off');
outputs = fsolve(@findpsialpha,guess,options,a,b,slen);
psi = outputs(1);
alpha = outputs(2);
[a1,b1,k,mode,a2,b2,k2] = solvebeambendingfin(psi,alpha,slen);
if abs(a2-a)<abs(a1-a)
    k = k2;
    mode = 1;
end



function errout = findpsialpha(in,ain,bin,slen)
tol = 1e-5;
psi = in(1);
alpha = in(2);
[a,b,k,mode, a2, b2,k2] = solvebeambendingfin(psi,alpha,slen);
if abs(ain-1)>tol
    errout1 = (ain - a)/(ain-1);
    errout2 = (ain - a2)/(ain - 1);
else
    errout1 = (ain - a)/tol;
    errout2 = (ain - a2)/tol;
end
errout3 = (bin - b)/bin;
errout4 = (bin - b2)/bin;
errtot1 = sqrt(errout1^2 + errout3^2);
errtot2 = sqrt(errout2^2 + errout4^2);
if errtot1 < errtot2
    errout(1) = errout1;
    errout(2) = errout3;
else
    errout(1) = errout2;
    errout(2) = errout4;
end

function [dx,dy,k,mode, a2, b2,k2] = solvebeambendingfin(psi,alpha,slen)
%Finds the value of k that matches psi and alpha.
%options = optimset('TolFun',1e-9,'MaxFunEvals',10000,'MaxIter',10000,'TolX',1e-15);%,'Display','off');
a2 = 0;
b2 = 0;
k2 = 0;
kmin = sqrt((cos(psi)+1)/2)*(1+eps);%.000000000000001;
if kmin == 0
    kmin = sqrt(eps);
end
if abs(sin(psi))<eps
    kmax = sqrt(1-eps);
    alphamax = findalpha(kmax,psi,1);
    alphamin = pi^2;
    if alpha >= alphamin & alpha <= alphamax
        k = fzero(@solvefnfin,[kmin kmax],[],psi,alpha,1);
        mode = 1;
    else
        k = 0;
        mode = 1;
    end
elseif sin(psi)>eps
    kmax = sqrt(1-eps);
    alphamax = findalpha(kmax,psi,1);
    alphamin = findalpha(kmin,psi,1);
    if alpha<=alphamax & alpha >= alphamin
        k = fzero(@solvefnfin,[kmin kmax],[],psi,alpha,1);
    else
        k=0;
    end
    mode = 1;
else
    alphalow1 = findalpha(kmin,psi,1);
    alphalow2 = findalpha(kmin,psi,2);
    if psi>0
        psi = psi - 2*pi;
    end
    kmax = fminbnd(@findalpha,kmin,1,[],psi,1);
    alphahigh1 = findalpha(kmax,psi,1);
    alphahigh2 = findalpha(sqrt(1-eps),psi,2);
    if alpha > alphahigh1
        k2 = fzero(@solvefnfin, [kmax sqrt(1-eps)],[],psi,alpha,1);
        [as2, bs2, s, theta] = guidednonlinsp(k2,psi*180/pi,1);
        sa = cos(psi)*(cos(theta)).^2 + sin(psi)*sin(theta).*cos(theta);
        sb = cos(psi)*cos(theta).*sin(theta) + sin(psi)*(sin(theta)).^2;
        da = slen*alpha.*trapz(s,sa);
        db = slen*alpha.*trapz(s,sb);
        a2 = as2 + da(end);
        b2 = bs2 + db(end);
    end
    if alpha>=min([alphalow1 alphahigh1]) && alpha <= max([alphalow1 alphahigh1])
        k = fzero(@solvefnfin,[kmin kmax],[],psi,alpha,1);
        mode = 1;
    elseif alpha >= min([alphalow2 alphahigh2]) && alpha <= max([alphalow2 alphahigh2])
        kmax = sqrt(1-eps);
        k = fzero(@solvefnfin,[kmin kmax],[],psi,alpha,2);
        mode = 2;
    else
        k = 0;
        mode = 1;
    end
end

%[alpha psi]
[as,bs,s,theta] = guidednonlinsp(k,psi*180/pi,mode);

if s==0
    dx = as;
    dy = bs;
else
    sa = cos(psi)*(cos(theta)).^2 + sin(psi)*sin(theta).*cos(theta);
    sb = cos(psi)*cos(theta).*sin(theta) + sin(psi)*(sin(theta)).^2;
    da = slen*alpha.*trapz(s,sa);
    db = slen*alpha.*trapz(s,sb);
    dx = as + da(end);
    dy = bs + db(end);
end

function errout = solvefnfin(k,psi,alpha,mode)
% Returns the error between alphaout, calculated using k and psi, and alpha
alphaout = findalpha(k,psi,mode);
errout = alphaout - alpha;
     

function alpha = findalpha(k,psi,mode)
if abs(sin(psi))<eps
    phi1 = 0;
elseif sin(psi)>eps
    phi1 = asin(sqrt((cos(psi) + 1)./(2*k.^2)));
else
    phi1 = -asin(sqrt((cos(psi) + 1)./(2*k.^2)));
end
if ~isreal(phi1)
    phi1 = real(phi1);
end
if mode==1
    phi2 = pi - phi1;
else
    phi2 = 2*pi + phi1;
end

f1 = elliptic12(phi1,k.^2);

f2 = elliptic12(phi2,k.^2);
alpha = (f2 - f1).^2;
        

function [a,b,s,theta] = guidednonlinsp(k,psi,mode)
%This function calculates end coordinates and loads for a beam.
%Inputs are k, the elliptic integral modulus, and psi, the angle of the
%applied force, and mode, an integer giving the bending mode. npts is the
%number of points to use in giving the beam shape and moment.
%Outputs are a, the nondimensional horizontal displacement, b, the
%nondimensional vertical displacment, and alpha, the nondimensional load.

npts = 101;
kmin = sqrt((cosd(psi) + 1)/2);
if k<=kmin
    a = NaN;
    b = NaN;
    theta = zeros(npts,1);
    s = theta;
    return
elseif k>=1
    a = NaN;
    b = NaN;
    theta = zeros(npts,1);
    s = theta;
    return
end
if abs(sind(psi))<eps
    phi1 = 0;
elseif sind(psi)>eps
    phi1 = asin(sqrt((cosd(psi) + 1)./(2*k.^2)));
else
    phi1 = -asin(sqrt((cosd(psi) + 1)./(2*k.^2)));
end
switch mode
    case 1
        phi2 = pi - phi1;%first mode
    case 2
        phi2 = 2*pi+phi1;%second mode
    case 3
        phi2 = 3*pi - phi1;%third mode
    case 4
        phi2 = 4*pi+phi1;%fourth mode
    case 5
        phi2 = 5*pi - phi1;%fifth mode
    case 6
        phi2 = 6*pi + phi1;%sixth mode
    case 7
        phi2 = 7*pi - phi1;%seventh mode
    case 8
        phi2 = 8*pi + phi1;%eigth mode
    case 9 %third mode, second type
        phi = asin(sqrt((cosd(psi) + 1)./(2*k.^2)));
        if sind(psi)>=0            
            phi1 = pi-phi;
            phi2 = 4*pi + phi;
        else
            phi1 = pi + phi;
            phi2 = 4*pi - phi;
        end
end
[f1,e1] = elliptic12(phi1,k.^2);
[f2,e2] = elliptic12(phi2,k.^2);
alpha = (f2 - f1).^2;
a = -(cosd(psi).*(f1 - f2 + 2*(e2 - e1)) + ...
    2*sind(psi).*k.*(cos(phi2) - cos(phi1)))./sqrt(alpha);
b = -(sind(psi).*(f1 - f2 + 2*(e2 - e1)) + ...
    2*k.*cosd(psi).*(cos(phi1)-cos(phi2)))./sqrt(alpha);
[s,theta] = guidedbeamshapesp(phi1,phi2,k,alpha,f1,npts);


function [s,theta] = guidedbeamshapesp(phi1,phi2,k,alpha,f1,npts)

if size(phi1,1)>1
    phi1 = phi1';
    phi2 = phi2';
    k = k';
    alpha = alpha';
    f1 = f1';
end

numpts = npts-1;
phi = zeros(numpts+1,length(phi1));
for counter = 1:length(phi1)
phi(:,counter) = (phi1(counter):(phi2(counter) - phi1(counter))/numpts:phi2(counter))';
end
temp2 = ones(size(phi,1),1);
[f2] = elliptic12(phi,temp2*(k.^2));
s = (f2 - temp2*f1)./sqrt(temp2*alpha);
theta = 2*(asin((temp2*k).*sin(phi))-asin((temp2*k).*sin(temp2*phi1)));

function [F,E,Z] = elliptic12(u,m,tol)
% ELLIPTIC12 evaluates the value of the Incomplete Elliptic Integrals 
% of the First, Second Kind and Jacobi's Zeta Function.
%
%   [F,E,Z] = ELLIPTIC12(U,M,TOL) where U is a phase in radians, 0<M<1 is 
%   the module and TOL is the tolerance (optional). Default value for 
%   the tolerance is eps = 2.220e-16.
%
%   ELLIPTIC12 uses the method of the Arithmetic-Geometric Mean 
%   and Descending Landen Transformation described in [1] Ch. 17.6,
%   to determine the value of the Incomplete Elliptic Integrals 
%   of the First, Second Kind and Jacobi's Zeta Function [1], [2].
%
%       F(phi,m) = int(1/sqrt(1-m*sin(t)^2), t=0..phi);
%       E(phi,m) = int(sqrt(1-m*sin(t)^2), t=0..phi);
%       Z(phi,m) = E(u,m) - E(m)/K(m)*F(phi,m).
%
%   Tables generating code ([1], pp. 613-621):
%       [phi,alpha] = meshgrid(0:5:90, 0:2:90);                  % modulus and phase in degrees
%       [F,E,Z] = elliptic12(pi/180*phi, sin(pi/180*alpha).^2);  % values of integrals
%
%   See also ELLIPKE, ELLIPJ, ELLIPTIC3, THETA, AGM.
%
%   References:
%   [1] M. Abramowitz and I.A. Stegun, "Handbook of Mathematical Functions", 
%       Dover Publications", 1965, Ch. 17.1 - 17.6 (by L.M. Milne-Thomson).
%   [2] D. F. Lawden, "Elliptic Functions and Applications"
%       Springer-Verlag, vol. 80, 1989

%   For support, please reply to 
%       moiseev[at]sissa.it, moiseev.igor[at]gmail.com
%       Moiseev Igor, 
%       34106, SISSA, via Beirut n. 2-4,  Trieste, Italy
%
%   The code is optimized for ordered inputs produced by the functions 
%   meshgrid, ndgrid. To obtain maximum performace (up to 30%) for singleton, 
%   1-dimensional and random arrays remark call of the function unique(.) 
%   and edit further code. 

if nargin<3, tol = eps; end
if nargin<2, error('Not enough input arguments.'); end

if ~isreal(u) || ~isreal(m)
    error('Input arguments must be real.')
end

if length(m)==1, m = m(ones(size(u))); end
if length(u)==1, u = u(ones(size(m))); end
if ~isequal(size(m),size(u)), error('U and M must be the same size.'); end

F = zeros(size(u)); 
E = F;              
Z = E;
m = m(:).';    % make a row vector
u = u(:).';

if any(m < 0) || any(m > 1), error('M must be in the range 0 <= M <= 1.'); end

I = uint32( find(m ~= 1 & m >eps) );
if ~isempty(I)
    [mu,J,K] = unique(m(I));   % extracts unique values from m
    K = uint32(K);
    mumax = length(mu);
    signU = sign(u(I));

    % pre-allocate space and augment if needed
	chunk = 7;
	a = zeros(chunk,mumax);
	c = a; 
	b = a;
	a(1,:) = ones(1,mumax);
	c(1,:) = sqrt(mu);
	b(1,:) = sqrt(1-mu);
	n = uint32( zeros(1,mumax) );
	i = 1;
	while any(abs(c(i,:)) > tol)                                    % Arithmetic-Geometric Mean of A, B and C
        i = i + 1;
        if i > size(a,1)
          a = [a; zeros(2,mumax)];
          b = [b; zeros(2,mumax)];
          c = [c; zeros(2,mumax)];
        end
        a(i,:) = 0.5 * (a(i-1,:) + b(i-1,:));
        b(i,:) = sqrt(a(i-1,:) .* b(i-1,:));
        c(i,:) = 0.5 * (a(i-1,:) - b(i-1,:));
        in = uint32( find((abs(c(i,:)) <= tol) & (abs(c(i-1,:)) > tol)) );
        if ~isempty(in)
          [mi,ni] = size(in);
          n(in) = ones(mi,ni)*(i-1);
        end
	end
     
    mmax = length(I);
	mn = double(max(n));
	phin = zeros(1,mmax);     C  = zeros(1,mmax);    
	Cp = C;  e  = uint32(C);  phin(:) = signU.*u(I);
	i = 0;   c2 = c.^2;
	while i < mn                                                    % Descending Landen Transformation 
        i = i + 1;
        in = uint32(find(n(K) > i));
        if ~isempty(in)     
            phin(in) = atan(b(i,K(in))./a(i,K(in)).*tan(phin(in))) + ...
                pi.*ceil(phin(in)/pi - 0.5) + phin(in);
            e(in) = 2.^(i-1) ;
            C(in) = C(in)  + double(e(in(1)))*c2(i,K(in));
            Cp(in)= Cp(in) + c(i+1,K(in)).*sin(phin(in));  
        end
	end
    
    Ff = phin ./ (a(mn,K).*double(e)*2);                                                      
    F(I) = Ff.*signU;                                               % Incomplete Ell. Int. of the First Kind
    Z(I) = Cp.*signU;                                               % Jacobi Zeta Function
    E(I) = (Cp + (1 - 1/2*C) .* Ff).*signU;                         % Incomplete Ell. Int. of the Second Kind
end

% Special cases: m == {0, 1}
m0 = find(m <= eps);
if ~isempty(m0), F(m0) = u(m0); E(m0) = u(m0); Z(m0) = 0; end

m1 = find(m == 1);
um1 = abs(u(m1)); 
if ~isempty(m1), 
    N = floor( (um1+pi/2)/pi );  
    M = find(um1 < pi/2);              
    
    F(m1(M)) = log(tan(pi/4 + u(m1(M))/2));   
    F(m1(um1 >= pi/2)) = Inf.*sign(u(m1(um1 >= pi/2)));
    
    E(m1) = ((-1).^N .* sin(um1) + 2*N).*sign(u(m1)); 
    
    Z(m1) = (-1).^N .* sin(u(m1));                      
end
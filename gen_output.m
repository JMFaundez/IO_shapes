clear all
close all

addpath('/scratch/josfa/Tools/matlab-tools/nek/')

mesh_input = 'mesh_11_BC.bc';

[nekdata,lr1,elmap,~,~,fields,emode,wdsz,etag,hdr] = readnek(mesh_input);

nelx = 200;
nely = 40;
[xx,yy,uu,vv] = reshapenek(nekdata,nelx,nely);

[ny, nx] = size(xx);

x0 = 0.14;
r_s = 0.01;

A = 0.2969;                                                                                           
B = 0.1260;                                                                                           
C = 0.3516;                                                                                           
D = 0.2843;                                                                                           
E = 0.1015;                                                                                           
                                                                                                     
naca = @(x) 5*0.08*(A*sqrt(x) - B*x - C*x.^2 + D*x.^3 - E*x.^4);                                      
dnaca = @(x) 5*0.08*(0.5*A*1./sqrt(x) - B - 2*C*x + 3*D*x.^2 - 4*E*x.^3); 

y0 = naca(x0);
%alpha = atan(dnaca(x_s));

shiox = uu*0;
shioy = vv*0;
shiox2 = uu*0;
shioy2 = vv*0;
[val, ind0] = min((xx(1,:)-x0).^2 + (yy(1,:)-y0).^2 );

xs = xx(1,ind0);
ys = yy(1,ind0);

Ntot = 0;
for i=1:nx
    d = sqrt((xx(1,i)-xs).^2 + (yy(1,i)-ys).^2);
    if d<=r_s
        alpha = atan(dnaca(xx(1,i)));
        cos_a = cos(alpha);
        sin_a = sin(alpha);
        dn = sqrt((xx(1,i)-xx(2,i)).^2 + (yy(1,i)-yy(2,i)).^2)
        shiox(1,i) = -cos_a/dn;
        shiox(2,i) = cos_a/dn;
        shioy(1,i) = -sin_a/dn;
        shioy(2,i) = sin_a/dn;
        shiox2(1,i) = 1;
        Ntot = Ntot + 1;
        if mod(i,lr1(1))==0
            Ntot = Ntot +1;
        end
    end
end
shiox = shiox/Ntot;
shioy = shioy/Ntot;
shiox2 = shiox2/Ntot;

meshdata = zeros(length(xx(:,1)),length(xx(1,:)),6);
meshdata(:,:,1) = xx;
meshdata(:,:,2) = yy;
meshdata(:,:,3) = shiox;
meshdata(:,:,4) = shioy;

nekdata_s = demeshnek(meshdata,lr1);

output_file = 'hw-y.fld';
writenek(output_file,nekdata_s,lr1,elmap,0,0,fields,emode,wdsz,etag);
alpha = atan(dnaca(xx(1,:)));
cos_a = cos(alpha);
sin_a = sin(alpha);
cos_a(yy(1,:)<0) = -cos_a(yy(1,:)<0);
for i=1:ny
    meshdata(i,:,3) = sin_a;
    meshdata(i,:,4) = cos_a;
end


nekdata_alpha = demeshnek(meshdata,lr1);

output_file = 'alpha_2d.fld';
writenek(output_file,nekdata_alpha,lr1,elmap,0,0,fields,emode,wdsz,etag);


meshdata(:,:,3) = shiox2;
meshdata(:,:,4) = shioy2;

nekdata_s2 = demeshnek(meshdata,lr1);

output_file = 'hw-y2.fld';
writenek(output_file,nekdata_s2,lr1,elmap,0,0,fields,emode,wdsz,etag);



ut0 = uu(1,:).*cos_a + vv(1,:).*sin_a;
ut1 = uu(2,:).*cos_a + vv(2,:).*sin_a;

%plot(xx(1,:),ut1-ut0)

% figure(1000)
% mesh(xx,yy,shiox)
% axis('equal')
% colorbar()
% view(2)

%status = writenek();

%figure(100)
%plot(x_s,y_s)
%axis('equal')
%

%figure(120)
%plot(xx(:,1),yy(:,1))
%
%
%X = reshape(xx,[ntot,1]);
%Y = reshape(yy,[ntot,1]);
%C = reshape(shiox,[ntot,1]);
%
%figure(200)
%hold on
%plot(x_s,y_s,'r-*')
%h = scatter(X,Y,3,C,'filled');
%
%%h.EdgeColor = 'none';
%axis('equal')





function [nekdata_s,dst0] = gen_2Dinput(x0,nekdata,SIGMAS,SIGMAN,nelx,nely,output_file)
%    clear all
%    close all
  	addpath matlab_script/
	
    %nelx = 200;
	%nely = 40;
	%ff_name = 'mesh_11_BC.bc';

	% Where?
	%x0 =  .15;
	N0 =  .00; % normalized by local displ. thickness
	%SIGMAS = 1.00; % normalized by local displ. thickness
	%SIGMAN = 0.50; % normalized by local displ. thickness

	% Output files
	%outfile = ['forcing/cdist_m',num2str(mesh_n),'_3d.fld'];

	%% Get GLL points and baseflow

	% read flow field
	%[nekdata,lr1,elmap,~,~,fields,emode,wdsz,etag,header] = readnek(ff_name);
	[nel,lr12,nfields] = size(nekdata);
	lr1 = sqrt(lr12);

	% reshape gll points
	[xx,yy,uu,vv,pp,TT] = reshapenek(nekdata,nelx,nely);

	clear nekdata

	% wall coordinates (by mesh construction)
	ds = sqrt( (xx(1,2:end)-xx(1,1:end-1)).^2 +...
			   (yy(1,2:end)-yy(1,1:end-1)).^2 );
	s  = [0 cumsum(ds)]; [~,is0] = min(xx(1,:)); s = s(is0) - s;

	nn = zeros(size(xx));
	ss = zeros(size(xx));
	nv = zeros([size(xx) 2]);

	% dx and dy distance from node to airfoil
	% nn normal distance, nv tangent and normal over the airfoil
	for i = 1:size(xx,2)
		dx = xx(:,i) - xx(1,i);
		dy = yy(:,i) - yy(1,i);
		
		ss(:,i) = s(i);
		nn(:,i) = sqrt(dx.^2 + dy.^2);
		
		n = [ dx(end); dy(end)]/nn(end,i);
		nv(:,i,1) = n(1); nv(:,i,2) = n(2);
	end

	% displacement thickness
	dst = zeros(size(xx,2),1);
 	U0 = 1;
 	pinf = 0;
 	pref = U0^2/2 + pinf;
 	Uv = sqrt((pref-pp(1,:))*2);
    ut = uu.*nv(:,:,2) - vv.*nv(:,:,1);

    for i=1:length(xx(1,:))
    	for j=1:length(xx(:,1))
    		if ut(j,i)>=0.99*abs(Uv(i))
    			integrand = 1 - ut(1:j,i)/Uv(i);
    			dst(i) = trapz(nn(1:j,i),integrand);
    			break
    		end
    	end
    end


%	for i = 1:size(xx,2)
%		nmax = max(nn(:,i)*0.8);
%		U = -nv(:,i,2) .* uu(:,i) + nv(:,i,1) .* vv(:,i);
%		
%		[Umax,jmax] = max(abs(U) .* (nn(:,i) < nmax));
%		jj = nn(:,i) > 1.5*nn(jmax,i);
%		%Uref = interp1(nn(jj,i),U(jj),nn(:,i),'spline');
%		Uref = U(end);
%		dst(i) = trapz(nn(:,i),(1 - U./Uref));
%	end

	%% Forcing field
	% disturbance location
	s0   = interp1(xx(1,:),s  ,x0);
	dst0 = interp1(xx(1,:),dst,x0);

	n0 = N0 * dst0;
	sigmas = SIGMAS * dst0;
	sigman = SIGMAN * dst0;
	%n0 = 0;
	%sigmas = SIGMAS;
	%sigman = SIGMAN;


	% forcing
	mask = exp(-(((ss-s0)/sigmas).^2 + ((nn-n0)/sigman).^2));
	%fs = (nn-n0)/dst0.*mask;
	%fn = (ss-s0)/dst0.*mask;
	fs = 0*mask;
	fn = mask;

	% global reference frame
	fx = nv(:,:,1).*fn - nv(:,:,2).*fs;
	fy = nv(:,:,2).*fn + nv(:,:,1).*fs;

	bx = fx;
	by = fy;

	meshdata = zeros(length(xx(:,1)),length(xx(1,:)),6);
	meshdata(:,:,1) = xx;
	meshdata(:,:,2) = yy;
	meshdata(:,:,3) = bx;
	meshdata(:,:,4) = by;
  
	nekdata_s = demeshnek(meshdata,lr1);

	%writenek(output_file,nekdata_s,lr1,elmap,0,0,fields,emode,wdsz,etag);
  
  %figure(100)
  %contourf(xx,yy,bx)
  %axis('equal')
  %colorbar()
end

function [nekdata_s,Ntot] = gen_2Doutput(x0,r_s,nekdata,nelx,nely, output_file)

	addpath matlab_script/
	
	%nelx = 300;
	%nely = 22;
	%mesh_input = 'fringe_m90.f00008';

	%[nekdata,lr1,elmap,~,~,fields,emode,wdsz,etag,hdr] = readnek(mesh_input);
	
	[nel,lr12,nfields] = size(nekdata);
	lr1 = sqrt(lr12);

	[xx,yy,uu,vv] = reshapenek(nekdata,nelx,nely);
	
	[ny, nx] = size(xx);
	
    A = 0.2969; 
	B = 0.1260;                                                                                         
	C = 0.3516;                                                                                         
	D = 0.2843;
	E = 0.1015;
	naca = @(x) 5*0.08*(A*sqrt(x) - B*x - C*x.^2 + D*x.^3 - E*x.^4);
	dnaca = @(x) 5*0.08*(0.5*A*1./sqrt(x) - B - 2*C*x + 3*D*x.^2 - 4*E*x.^3); 
	
	y0 = naca(x0);
	
	shiox = uu*0;
	shioy = vv*0;
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
	        dn = sqrt((xx(1,i)-xx(2,i)).^2 + (yy(1,i)-yy(2,i)).^2);
	        shiox(1,i) = -cos_a/dn;
	        shiox(2,i) = cos_a/dn;
	        shioy(1,i) = -sin_a/dn;
	        shioy(2,i) = sin_a/dn;
	        Ntot = Ntot + 1;
	        if mod(i,lr1)==0
	            Ntot = Ntot +1;
	        end
	    end
	end
%	shiox = shiox/Ntot;
%	shioy = shioy/Ntot;
	meshdata = zeros(length(xx(:,1)),length(xx(1,:)),6);
	meshdata(:,:,1) = xx;
	meshdata(:,:,2) = yy;
	meshdata(:,:,3) = shiox;
	meshdata(:,:,4) = shioy;
	
	nekdata_s = demeshnek(meshdata,lr1);

%	writenek(output_file,nekdata_s,lr1,elmap,0,0,fields,emode,wdsz,etag);
	
end

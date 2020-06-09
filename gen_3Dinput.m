function gen_3Dinput(x0,z0,SIGMAS, SIGMAN,SIGMAZ,output_file)
	addpath matlab_script/
	
	input3d = 'mesh_9_3D_fringe.bc';
	input2d = 'fringe_m90.f00008';
	nelx = 300;
	nely = 22;
	nelz = 14;
	
	%x0 = 0.15;
	%z0 = 0;
	%rs = 0.005;
	
	
	[threshold,lr2,elmap,~,~,fields,emode,wdsz,etag,hdr] = readnek(input2d);
	nxy = nelx*nely;
	
	[data3d,lr1,elmap,time,istep,fields,emode,wdsz,etag,header,status] = readnek(input3d);
	
	data_shio = data3d;
	
	gll2 = lr1(1)^2;
	Ntot = 0;
	for k=1:nelz
		ik = linspace((k-1)*nxy+1,k*nxy,nxy);
		for gl=1:lr1(1)
			nodes = linspace((gl-1)*gll2+1,gl*gll2,gll2);
			zi = max(max(data_shio(ik,nodes,3)));
			dz = abs(zi-z0);
	        [nekdata2,dst0] = gen_2Dinput(x0,threshold,SIGMAS,SIGMAN,nelx,nely,'...');
	        sigmaz = SIGMAZ*dst0;
	  		data_shio(ik,nodes,4) = nekdata2(:,:,3)*exp(-(dz/sigmaz)^2);  
	  		data_shio(ik,nodes,5) = nekdata2(:,:,4)*exp(-(dz/sigmaz)^2);
	  		data_shio(ik,nodes,6) = 0;  
	  		data_shio(ik,nodes,8) = threshold(:,:,4);  
		end
	end
	
	%output_file = 'hw-y1.fld';
	status = writenek(output_file,data_shio,lr1,elmap,time,istep,fields,emode,wdsz,etag)
end

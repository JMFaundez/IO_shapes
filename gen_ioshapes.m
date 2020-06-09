clear all 
close all


% gen ouput shapes
x0 = 0.1;
z0 = [-0.005, 0.005];
rs= 0.002;

for i=1:length(z0)
%    filename = ['hw-y',num2str(i),'.fld'];
%    gen_3Doutput(x0,z0(i),rs,filename)
end


% gen input shapes
x0 = 0.2;
z0 = 0;
SIGMAS = 3;
SIGMAN = 5;
SIGMAZ = 1.5;

for i=1:length(z0)
    filename = ['plasma',num2str(i),'.fld'];
    gen_3Dinput(x0,z0,SIGMAS, SIGMAN,SIGMAZ,filename)
end

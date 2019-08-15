function [K,ScanParams]=ReadFID(fn)

[Enc,AntiAlias,Direction,scaling_read]=ReadMethods(fn);

fid = fopen([fn,filesep,'fid']);
R = fread(fid,'int32','l');
fclose(fid);
 
K1=reshape(R,[Enc(1)*2,Enc(2),Enc(3)]);
K1=permute(K1,[2 3 1]);
for i=1:Enc(2); for j=1:Enc(3); for k=1:2:Enc(1)*2
KC1(i,j,ceil(k/2))=K1(i,j,k)+1i*K1(i,j,k+1);            
        end
    end
end

if scaling_read==-1
disp('Reversing Readout Direction');
KC1=KC1(:,:,end:-1:1); 
KC1=circshift(KC1,1,3); %31-7-2018 frequency shift 
end


% If FH 
if strcmp(Direction,'H_F')==1
    KC1=permute(KC1,[3 2 1]); 
end

K=permute(KC1,[3 1 2]); %%???

ScanParams.Enc=Enc; 
ScanParams.AntiAlias=AntiAlias; 
ScanParams.Direction=Direction;
ScanParams.scaling_read=scaling_read; 
end
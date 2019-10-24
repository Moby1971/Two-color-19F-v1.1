function [K,ScanParams]=ReadFIDSequentialCS(fn,ndirs,dirs,ky,kz)
[Enc,AntiAlias,Direction,scaling_read]=ReadMethods(fn);

fid = fopen([fn,filesep,'fid']);
R = fread(fid,'int32','l');
fclose(fid);
% to do: add sequential option here... 
% Should also take directions into account -or add these before in options
%  P.Seq.Sequential=1
%  P.Seq.Directions={'FH','HF','LR','RL'}

 
K1=reshape(R,[2,Enc(1),ndirs,Enc(2),Enc(3)]); %[2,ndirs,kx,ky,kz]   DATA BEFORE SORTING
K1=permute(K1,[1 3 4 5 2]);

K1=reshape(K1,[2,ndirs,Enc(2)*Enc(3),Enc(1)]);    % [2, ndirs, ky*kz, kx]

KC1 = zeros(Enc(2),Enc(3),Enc(1),ndirs);   % [ky, kz, kx, ndirs]

% SORT THE DATA
% KSPACE TRAJECTORY IS DEFINED IN    ky and kz

for w = 1 : Enc(2)*Enc(3)   % ky and kz directions
 
        for k=1:Enc(1) % readout = kx
            for q=1:ndirs
                KC1(ky(w),kz(w),k,q) = KC1(ky(w),kz(w),k,q) + K1(1,q,w,k)+1i*K1(2,q,w,k);
            end
        end
        
end



for jj=1:ndirs
    Ktemp=KC1(:,:,:,jj);
    Direction=dirs{jj};
    if strcmp(Direction,'LR') || strcmp(Direction,'HF')
        scaling_read=1;
    else
        scaling_read=-1;
    end
    
    if scaling_read==-1
        Ktemp=Ktemp(:,:,end:-1:1);
        Ktemp=circshift(Ktemp,1,3); %31-7-2018 frequency shift
    end
    
    % If FH
    if strcmp(Direction,'HF')==1 || strcmp(Direction,'FH')==1
        Ktemp=permute(Ktemp,[3 2 1]);
    end

    K{jj}=permute(Ktemp,[3 1 2]); %%???
    
    ScanParams{jj}.Enc=Enc;
    ScanParams{jj}.AntiAlias=AntiAlias;
    ScanParams{jj}.Direction=Direction;
    ScanParams{jj}.scaling_read=scaling_read;
end

end
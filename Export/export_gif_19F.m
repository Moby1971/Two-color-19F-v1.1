function export_gif_19F(gifexportpath,tag,PFCEmap,PFOBmap,PFCEcmap,PFOBcmap,image_orient)

% Exports PFCE and PFOB maps to animated gif


PFOBmap = permute(PFOBmap,image_orient);
PFCEmap = permute(PFCEmap,image_orient);

[number_of_images,dimx,dimy] = size(PFCEmap);

% increase the size of the matrix to make the exported images bigger

numrows = 2*dimx;
numcols = 2*dimy;

delay_time = 5/number_of_images;  % show all gifs in 5 seconds



% Export the PFOB maps to gifs

for idx = 1:number_of_images
    
    % because the color maps are arrays of size 64
    % the t2map needs to be mapped onto the range of [0, 63] and cast in an
    % unsigned integer 8 for gif export
    image = uint8(round((63/max(PFOBmap(:)))*resizem(squeeze(PFOBmap(idx,:,:)),[numrows numcols])));
    
    if idx == 1
        imwrite(image,PFOBcmap,[gifexportpath,'/PFOB-',tag,'.gif'],'DelayTime',delay_time,'LoopCount',inf);
    else
        imwrite(image,PFOBcmap,[gifexportpath,'/PFOB-',tag,'.gif'],'WriteMode','append','DelayTime',delay_time);
    end
    
end


% Export the PFCE maps to GIF

for idx = 1:number_of_images
    
    % because the color maps are arrays of size 64
    % the t2map needs to be mapped onto the range of [0, 63] and cast in an
    % unsigned integer 8 for gif export
    image = uint8(round((63/max(PFCEmap(:)))*resizem(squeeze(PFCEmap(idx,:,:)),[numrows numcols])));
    
    if idx == 1
        imwrite(image,PFCEcmap,[gifexportpath,'/PFCE-',tag,'.gif'],'DelayTime',delay_time,'LoopCount',inf);
    else
        imwrite(image,PFCEcmap,[gifexportpath,'/PFCE-',tag,'.gif'],'WriteMode','append','DelayTime',delay_time);
    end
    
end


end
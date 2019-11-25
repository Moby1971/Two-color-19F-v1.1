function export_gif_19F(app,gifexportpath,tag,C1map,C2map,C1cmap,C2cmap)

% Exports 19F iamges to animated gif


[number_of_images,dimx,dimy] = size(C1map);

% increase the size of the matrix to make the exported images bigger

numrows = 4*dimx;
numcols = 4*dimy;

delay_time = 5/number_of_images;  % show all gifs in 5 seconds



% Export the C2 maps to gifs

for idx = 1:number_of_images
    
    % because the color maps are arrays of size 64
    % the t2map needs to be mapped onto the range of [0, 63] and cast in an
    % unsigned integer 8 for gif export
    image = rot90(uint8(round((63/max(C2map(:)))*resizem(squeeze(C2map(idx,:,:)),[numrows numcols]))));
    
    if idx == 1
        imwrite(image,C2cmap,[gifexportpath,filesep,app.c2name,'-',tag,'.gif'],'DelayTime',delay_time,'LoopCount',inf);
    else
        imwrite(image,C2cmap,[gifexportpath,filesep,app.c2name,'-',tag,'.gif'],'WriteMode','append','DelayTime',delay_time);
    end
    
end


% Export the C1 maps to GIF

for idx = 1:number_of_images
    
    % because the color maps are arrays of size 64
    % the t2map needs to be mapped onto the range of [0, 63] and cast in an
    % unsigned integer 8 for gif export
    image = rot90(uint8(round((63/max(C1map(:)))*resizem(squeeze(C1map(idx,:,:)),[numrows numcols]))));
    
    if idx == 1
        imwrite(image,C1cmap,[gifexportpath,filesep,app.c1name,'-',tag,'.gif'],'DelayTime',delay_time,'LoopCount',inf);
    else
        imwrite(image,C1cmap,[gifexportpath,filesep,app.c1name,'-',tag,'.gif'],'WriteMode','append','DelayTime',delay_time);
    end
    
end


end
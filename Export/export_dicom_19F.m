function export_dicom_19F(directory,im,parameters,tag)


folder_name = [directory,[filesep,'DICOM-',tag]];
if (~exist(folder_name, 'dir')); mkdir(folder_name); end
delete([folder_name,filesep,'*']);

[nr_images,dimx,dimy] = size(im);

dcmid = dicomuid;   % unique identifier
dcmid = dcmid(1:50);

for i=1:nr_images
    dcm_header = generate_dicomheader_19F(parameters,i,dimx,dimy,dcmid);
    fn = ['0000',num2str(i)];
    fn = fn(size(fn,2)-4:size(fn,2));
    fname = [directory,filesep,'DICOM-',tag,filesep,tag,'_',fn,'.dcm'];
    image = rot90(squeeze(cast(round(im(i,:,:)),'uint16')));
    dicomwrite(image, fname, dcm_header);
end


end
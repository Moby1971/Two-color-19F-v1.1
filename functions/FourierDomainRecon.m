function R= FourierDomainRecon(params)

k1=params.k1;
k2=params.k2; 
k3=params.k3;
k4=params.k4; 
nx=params.nx;
ny=params.ny; 
A1=params.A1; 
A2=params.A2;
A3=params.A3;
A4=params.A4; 
FS=params.FS;
F=params.F;
N=params.N; 
rr1 = @(I) reshape(I,[nx,ny]);

%making operators...
Z=zeros(nx,ny); Z(1,1)=1; 
ConvUnitResponse=F*Z(:); 


if params.discreteFTresponse
    ConvResponse1=A1*Z(:);
    ConvResponse2=A4*Z(:);
    ConvResponse3=A3*Z(:);
    ConvResponse4=A2*Z(:);
    
    FConvResponse1=F*(ConvResponse1);
    FConvResponse2=F*(ConvResponse2);
    FConvResponse3=F*(ConvResponse3);
    FConvResponse4=F*(ConvResponse4);
    
else
    %% TEMP TEMP TEMP
    BW=4.4643e+04; %(scan 13& 14)
    BWpix=BW/128;
    ppm=282.5685;
    [PFCE,PFCE_alpha,PFOB,PFOB_alpha]=calcspectra_BW(ppm,BWpix); %126.4??
    PFOB_alpha=[345 503 340 900 340];
    PFOB_alpha=PFOB_alpha./sum(PFOB_alpha(:));
    
    omega=linspace(0,1-(1/nx),nx); 
    for ii=1:length(PFOB_alpha) 
        response(ii,:)=(ones(1,nx).*PFOB_alpha(ii)).*exp(-1i*2*pi*omega*(-(PFOB(ii))));
    end
    response_halfshift=(ones(1,nx)).*exp(-1i*2*pi*omega*(0));

    response1D=sum(response,1);
    response2D=ones(nx,ny)./sqrt(nx*ny);
%     response2D=bsxfun(@times,response2D,response_halfshift);

    
FConvResponse1=bsxfun(@times,response2D,response1D.');
FConvResponse2=bsxfun(@times,response2D,conj(response1D));
FConvResponse3=bsxfun(@times,response2D,response1D);
FConvResponse4=bsxfun(@times,response2D,response1D');

% ConvUnitResponse=F*Z(:);  %without shift (for now)
ConvUnitResponse=ones(nx,ny)./sqrt(nx*ny);

%% temp comparison 
    discreteresponse3=(rr1(F*A3*Z(:)));
    discreteresponse4=(rr1(F*A2*Z(:)));

    
    figure(101);clf; imshow(abs(discreteresponse3),[]);
    
    figure(102);clf;  hold on; plot(abs(discreteresponse3(1,:)));
    plot(abs(FConvResponse3(1,:)));hold off;  title('abs response 3')
    
        figure(103); clf;  hold on; plot(angle(discreteresponse3(1,:)));
    plot(angle(FConvResponse3(1,:)));hold off; title('angle response 3')
    
        figure(104); clf;  hold on; plot(angle(discreteresponse3(:,1)));
    plot(angle(FConvResponse3(:,1)));hold off; title('angle response 3 - other dir')

    figure(105);clf;  hold on; plot(abs(discreteresponse4(:,1)));
    plot(abs(FConvResponse4(:,1)));hold off; title('abs response 4')


end


%% 
opUnitResponse=opDiag(ConvUnitResponse(:));
opFConvResponse1=opDiag(FConvResponse1(:));
opFConvResponse2=opDiag(FConvResponse2(:));
opFConvResponse3=opDiag(FConvResponse3(:));
opFConvResponse4=opDiag(FConvResponse4(:));

switch params.nrmeas
    case 2
        bigC=[opUnitResponse opFConvResponse3; opUnitResponse opFConvResponse4];
        measurements=[k3(:);k4(:)];
    case 4
        bigC=[opUnitResponse opFConvResponse1; opUnitResponse opFConvResponse2;...
        opUnitResponse opFConvResponse3;opUnitResponse opFConvResponse4];
         measurements=[k1(:);k2(:);k3(:); k4(:)];
end

% inverting matrix..
Ks=bigC\measurements;

lambda=1e4
niter=20
xCG=nl_conjgrad_fluor_kspace(bigC, measurements, Ks,niter,Ks,lambda,nx,ny*2,1,FS)


Ksolve1=squeeze(Ks(1:N)); 
Ksolve2=squeeze(Ks(N+1:2*N)); 
Ksolve1(isnan(Ksolve1))=0;
Ksolve2(isnan(Ksolve2))=0;

KsolveCG1=squeeze(xCG(1:N)); 
KsolveCG2=squeeze(xCG(N+1:2*N)); 
KsolveCG1(isnan(KsolveCG1))=0;
KsolveCG2(isnan(KsolveCG2))=0;

R=[Ksolve1,Ksolve2];

figure(10);
subplot(221); 
imshow(abs(rr1(FS'*Ksolve1)),[]); title('non-iterative fourier reconstruction')
subplot(222); 
imshow(abs(rr1(FS'*Ksolve2)),[])

subplot(223); 
imshow(abs(rr1(FS'*KsolveCG1)),[]); title('CG fourier reconstruction')
subplot(224); 
imshow(abs(rr1(FS'*KsolveCG2)),[])
end
% The results of CaImAn are saved in an estimates object. This is stored inside the cnmf object, i.e. it can be accessed using cnmf.estimates. The variables of interest are:
% 
% estimates.A: Set of spatial components. Saved as a sparse column format matrix with dimensions (# of pixels X # of components). Each column corresponds to a spatial component.
% 
% 
% estimates.b: Set of background spatial components (for 2p analysis): Saved as a numpy array with dimensions (# of pixels X # of components). Each column corresponds to a spatial background component.
% 
% estimates.f: Set of temporal background components (for 2p analysis). Saved as a numpy array with dimensions (# of background components X # of timesteps). Each row corresponds to a temporal background component.
% 
% estimates.S: Deconvolved neural activity (spikes) for each component. Saved as a numpy array with dimensions (# of background components X # of timesteps). Each row corresponds to the deconvolved neural activity for the corresponding component.
% 
% estimates.YrA: Set of residual components. Saved as a numpy array with dimensions (# of components X # of timesteps). Each row corresponds to the residual signal after denoising the corresponding component in estimates.C.
% 


FileName='Slide1.sld - 1-pwer9-Frame2-noCNO-reversespm-091922_E_10_Iter_3500_output.hdf5';
PathName='/Users/restrepd/Documents/Projects/SFTP/20220919-mmHippoCamKGrin1/file1 for raw traces check/';

figureNo=0;

% estimates.b: Set of background spatial components (for 2p analysis): Saved as a numpy array with dimensions (# of pixels X # of components). Each column corresponds to a spatial background component
b=h5read([PathName FileName],'/estimates/b' );

% estimates.f: Set of temporal background components (for 2p analysis). Saved as a numpy array with dimensions (# of background components X # of timesteps). Each row corresponds to a temporal background component
f=h5read([PathName FileName],'/estimates/f' );
figureNo=figureNo+1;
figure(figureNo)
plot(f(:,1))



% estimates.YrA: Set of residual components. Saved as a numpy array with dimensions (# of components X # of timesteps). Each row corresponds to the residual signal after denoising the corresponding component in estimates.C.
% YrA=h5read([PathName FileName],'/estimates/YrA' );


% estimates.C: Set of temporal components. Saved as a numpy array with dimensions (# of components X # of timesteps) 
% Each row corresponds to a temporal component denoised and deconvolved.
C=h5read([PathName FileName],'/estimates/C' );
figureNo=figureNo+1;
figure(figureNo)
plot(C(:,1))

% estimates.F_dff: Set of DF/F normalized temporal components. Saved as a numpy array with dimensions (# of components X # of timesteps). 
% Each row corresponds to the DF/F fluorescence for the corresponding component
dFF=h5read([PathName FileName],'/estimates/F_dff' );
figureNo=figureNo+1;
figure(figureNo)
plot(dFF(:,1))

% estimates.S: Deconvolved neural activity (spikes) for each component. 
% Saved as a numpy array with dimensions (# of background components X # of timesteps). Each row corresponds to the deconvolved neural activity for the corresponding component
S=h5read([PathName FileName],'/estimates/S' );
figureNo=figureNo+1;
figure(figureNo)
plot(S(:,1))

pfft=1;


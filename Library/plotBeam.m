

function plotBeam( data )

% N  =  size(data,2) / 4;
% 
% subplot(2,2,1);  imagesc( data(:,0*N+1 : 1*N)); axis image off; colorbar; title('Porosity');
% subplot(2,2,2);  imagesc( data(:,1*N+1 : 2*N)); ); axis image off; colorbar; title('c [m/s]');
% subplot(2,2,3);  imagesc( data(:,2*N+1 : 3*N)); ); axis image off; colorbar; title('alpha [dB/(MHz*cm)]');
% subplot(2,2,4);  imagesc( data(:,3*N+1 : 4*N)); ); axis image off; colorbar; title('Pressure [a.u.]');

bar(data);







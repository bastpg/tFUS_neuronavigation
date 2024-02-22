function [projections]=project(vS,fS,vT,fT)

% id = any(isnan(vS),2);
% fS(id,:) = [];
% id2 = any(isnan(vT),2);
% fT(id,:) = [];
% vS(any(isnan(vS),2),:) = [];
% vT(any(isnan(vT),2),:) = [];

TRS = triangulation(fS,vS); 
normalsS=vertexNormal(TRS);


[IDXsource,Dsource]=knnsearch(vT,vS);
vector_s_to_t=vT(IDXsource,:)-vS;

projections=vS+[(sum(vector_s_to_t.*normalsS,2)./(norm(normalsS).^2)).*normalsS(:,1) (sum(vector_s_to_t.*normalsS,2)./(norm(normalsS).^2)).*normalsS(:,2) (sum(vector_s_to_t.*normalsS,2)./(norm(normalsS).^2)).*normalsS(:,3)];

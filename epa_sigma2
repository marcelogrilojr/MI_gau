% function [H,S,D]=mi_epa(X,n,p)
function [H]=mi_epa(X,n,p)

 u=X;
 [canais,N]=size(u);

%media nula para U
 u(:,:)=u(:,:)-kron(mean(u(:,:)')',ones(1,N*1)); 
 [UU,S,VV]=svd(u(:,:)',0);
 Q= pinv(S)*VV';

%normalizacao
for c=1:canais
   z(c,:)=u(c,:)/max(u(c,:)); 
end


%%               
%Epanechnikov 
%% PARTE LINEAR
%Inicializacao da matriz B com valores aleatorios
    B=eye(canais);
%   Bini=B;
    I=eye(canais);
        
    %Para 19 canais
    
    sigma_epa = [5.3 4.9 5 4.7 5.2 5 4.9 5.1 5.3 5.1 5.3 5.2 4.9 5.1 5 5.1 4.9 5.2 5.2 5.1 4.8 4.9 5.3 5.1 4.9 5 5.3 4.7 5.2 5 4.8 5.1]; 
    
    %Passo de adaptacao linear 
    mu=0.001;
        
    %tamanho da janela
    a=100;
    b=a-1;
%Janela Deslizante
    for k=1:N-b
    
        y=B*z;
            %vetor y da janela deslizante    
            yj=y(:,k:k+b);
            zj=z(:,k:k+b);
                     
            q=length(yj);
            dH_yb = zeros(canais,1,canais);
                            
            %kernel
            epa_yj = zeros(canais,a);
%             epa_zj = zeros(canais,a);
                                        
            for r = 1:canais
            %Estimativa da pdf com kernel Epanechnikov
                for i=1:q
%                     cont_Y=0;
%                     cont_Z=0;
                    for j=1:q
                        
                        %Condicao para kernel Epanechnikov                            
                        ysigma=yj(r,i)-yj(r,j);                     
                        if (abs(ysigma)<=sigma_epa(r))
                            %epa_yj(r,i)=epa_yj(r,i)+(1/(sqrt(2*pi)*sigma_epa(1,r)))*exp((-((yj(r,i)-yj(r,j))^2))/(2*sigma_epa(1,r)^2));
                            epa_yj(r,i)=epa_yj(r,i)+((3/(4*sigma_epa(r)))*(1-((ysigma/sigma_epa(r))^2)));                                   
                            %cont_Y=cont_Y+1;
                        end
                        
%                         zsigma=zj(r,i)-zj(r,j);
%                         if (abs(zsigma)<=sigma_epa(r))
%                             %epa_zj(r,i)=epa_zj(r,i)+(1/(sqrt(2*pi)*sigma_epa(1,r)))*exp((-((zj(r,i)-zj(r,j))^2))/(2*sigma_epa(1,r)^2));
%                             epa_zj(r,i)=epa_zj(r,i)+((3/(4*sigma_epa(1,r)))*(1-((zsigma/sigma_epa(1,r))^2)));
%                             cont_Z=cont_Z+1;
%                         end
                        
                    end
                    epa_yj(r,i)=epa_yj(r,i)./q;
%                     epa_zj(r,i)=epa_zj(r,i)./cont_Z;
                end
            end
                        
            for r = 1:canais
                for i=1:q
                    for j=1:q
                        
                    %Condicao para kernel Epanechnikov                            
                    ysigma=yj(r,i)-yj(r,j);  
                    
                        if (abs(ysigma)<=sigma_epa(r))
                            %Derivada da Entropia H(Y,B)
                            dH_yb(:,:,r)=dH_yb(:,:,r)+((1/q^2)*(1/epa_yj(r,i))*(6/(4*sigma_epa(1,r)))*((ysigma)/(sigma_epa(1,r)^2))*(zj(:,i)-zj(:,j)));
                        end
                                   
                    end%for
                end
            end

    %SF=zeros(canais,canais);
    for t=1:canais
        SF(:,t)= dH_yb(:,:,t)';
    end
    
    % Gradiente Relativo:
    B =(I-mu*SF)*B;
        
     for c=1:canais
        B(c,:)=B(c,:)./norm(B(c,:));
     end

    end
%%
    %Recuperacao do Sinal
     S=B*z;
  
% %%
%     %PARA O EEGLAB
%     % Matriz de mistura no SOBI
%    H = pinv(Q)*B; 
%    % H=B
%     [vcp,D] = eig(real(B*B'));
%     S=y;
% %     S=B'*z(:,:) % estimated source activities
%     %save('PNL_SOBI');
% % end

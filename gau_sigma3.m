function [S]=mi_gau(X,n,p)
% function [H,S,D]=mi_gau(X,n,p)
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
%  z=Z;
%x=u;

%%               
%Epanechnikov 
%Inicializacao da matriz B com valores aleatorios
    B=eye(canais);
%   Bini=B;
    I=eye(canais);
        
    %Para 19 canais
    sigma = [7.3 7.5 7.1 7.4 7.2 7.5 7.2 7.1 7.3 7.1 7.3 7.5 7.4 7.1 7.5 7.1 7.3 7.2 7.1 7.4 7 7.1 7.4 7.2 7.1 7.3 7.4 7 7.2 7.1 7.2 7.3]; 
        
    mu=0.01;
        
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
            gau_yj_NL_c = zeros(canais,a);
%             gau_zj_NL_c = zeros(canais,a);
                  
            for r = 1:canais
            %Estimativa da pdf com kernel gaussiano
                for i=1:q
                    cont_Y=0;
%                     cont_Z=0;
                    
                    for j=1:q
                        gau_yj_NL_c(r,i)=gau_yj_NL_c(r,i)+(1/(sqrt(2*pi)*sigma(1,r)))*exp((-((yj(r,i)-yj(r,j))^2))/(2*sigma(1,r)^2));
                        cont_Y=cont_Y+1;

%                         gau_zj_NL_c(r,i)=gau_zj_NL_c(r,i)+(1/(sqrt(2*pi)*sigma(1,r)))*exp((-((zj(r,i)-zj(r,j))^2))/(2*sigma(1,r)^2));
%                         cont_Z=cont_Z+1;
                    end
                    gau_yj_NL_c(r,i)=gau_yj_NL_c(r,i)./cont_Y;
%                     gau_zj_NL_c(r,i)=gau_zj_NL_c(r,i)./cont_Z;
                end
            end
                        
            for r = 1:canais
                for i=1:q
                    for j=1:q
                    %Derivada da Entropia H(Y,B)---LINEAR
                    %dH_yb(:,:,r)=dH_yb(:,:,r)+((1/q^2)*(1/gau_yj_NL_c(r,i))*(1/(sqrt(2*pi)*sigma(1,r)))*(1/sigma(1,r)^2)*exp((-((yj(r,i)-yj(r,j))^2))/(2*sigma(1,r)^2))*(-(yj(r,i)-yj(r,j)))*(zj(:,i)-zj(:,j)));
                    dH_yb(:,:,r)=dH_yb(:,:,r)+((1/q^2)*(1/gau_yj_NL_c(r,i))*(1/(sqrt(2*pi)*sigma(1,r)))*(1/sigma(1,r)^2)*exp((-((yj(r,i)-yj(r,j))^2))/(2*sigma(1,r)^2))*((yj(r,i)-yj(r,j)))*(zj(:,i)-zj(:,j)));

                    end
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
%     y_epa=y

% %%
%     %PARA O EEGLAB
%     % Matriz de mistura no SOBI
%     H = pinv(Q)*B; 
%    % H=B
%     [vcp,D] = eig(real(B*B'));
%     S=y;
%     %save('PNL_SOBI');
% % end

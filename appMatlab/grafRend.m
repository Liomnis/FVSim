function [] = grafRend(G_tr,T_tr,Rs,Rp,Vocn,Iscn,Kv,Ki,...
        Impp,Ns,Np,mg,largoPV,anchoPV,modelo)
%% entrada de valores de Temp e Irradiancia
T = T_tr + 273.15; % temp. de trabajo
Tn = 298.15; % temp. de referencia 25+273=298
dT = T-Tn;
G = G_tr; % irradiancia de trabajo del panel
Gn = 1000; % irradiancia de referencia 1000 W/m2
%% constantes físicas
q = 1.60e-19 ; % carga del Electrón
K = 1.38e-23 ; % constante de Boltzman
%% condiciones iniciales
p = 2.2;
Vt = (Ns*K*T)/q;
Io = (Iscn+Ki*dT)/(exp((Vocn+Kv*dT)/Vt)-1);% Io1 = Io2
%% ecuaciones preliminares
tT = 1;
Ipv = (Iscn+Ki*dT)*(G/Gn)*tT;
cont = 1;
%%
paso = .1;
h = waitbar(0,'Construyendo gráficas, por favor espere...');
for V = 0: paso: Vocn
    waitbar(V/Vocn) %waitbar
    x = Impp; % valor inicial de I
    Tol = 0.001; % tolerancia o error
    count = 0; % contador de iteraciones
    error = 1; %derivada de la ftn
    %% ecuaciones preliminares
%     f1 = Ipv-Io*(exp((V+I*Rs)/Vt)+exp((V+I*Rs)/((p-1)*Vt))+2)-((V+I*Rs)/Rp)-I;% Ipv(eq3)
        f1 = @(x) Ipv-Io*(exp((V+x*Rs)/Vt)+exp((V+x*Rs)/((p-1)*Vt))+2)-...
            ((V+x*Rs)/Rp)-x;% Ipv(eq3)
%     fini = inline(f1);% convierto ftn a inline
    f =  f1(x);% evaluo f para el valor inicial de I
    fini_deriv = diff(sym(f1),1);
    %%
    fprintf('iter       V           I          f(I)       error \n')
    while  (error > Tol || abs(f) > Tol)
        count = count + 1;
        fprime = eval(fini_deriv);
        xnew = x - (f/fprime);%
        error = abs(x-xnew);%
        x = xnew;
        f = f1(x);% evaluo el nuevo valor de f(I)
        fprintf('%3i%12.4f%12.4f%12.4f%12.4f \n',count,V,x,f,error)
    end
    I_array(1,cont) = x;
    V_array(1,cont) = V;
    P = V * x;
    P_array(1,cont) = P;
    cont = cont + 1;
    %
    if x < 0
        break
    end
end
close(h) % closed the waitbar

%% guardar valores de salida V I P en fichero txt
GT = ({[num2str(round(T_tr)),'_',num2str(round(G_tr))]});
GT = char(GT);
results_PIV = [P_array;I_array;I_array];
% modelostr = char(modelo);
% filetxt = char({['PIV_',modelo,'.txt']});
filetxt = char({['GrafRend_PIV_',modelo,'_',GT,'.txt']});
fileID = fopen(filetxt,'w');
fprintf(fileID,'%8s %8s %8s\n','P(W)', 'V(V)', 'I(A)');
fprintf(fileID,'%8.2f %8.2f %8.2f\n',results_PIV);
fclose(fileID);

%%
p=P_array;
p = find(p==max(P_array));
tmpp = V_array(p);
cmpp = I_array(p);

    %% grafico
    figure(mg+1)
    %% grafica de relleno P ideal
%     xi=linspace(0,max(V_array),max(V_array));
%     fi=@(xi) 0*xi+max(I_array);
%     yi=fi(xi);
%     xi0=0;
%     xi1=max(V_array);
%     xix=[xi0 xi0 xi(xi>xi0 & xi<xi1) xi1 xi1];
%     yiy=[0 fi(xi0) yi(xi>xi0 & xi<xi1) fi(xi1) 0];
%     fill(xix,yiy,'c');
%     hold on
    %%
    plot(V_array,I_array,'r');
%     set(gca,'xtick', 0:2:max(V_array));% valores del eje X
%     set(gca,'ytick', [-3 -1 2 4 6 10 15]);
    axis([0 max(V_array)+1 0 max(I_array)+1]);
    xlabel('Tensión (V)');
    ylabel('Corriente (A)');
    title('Rendimiento del panel');
    text(1,max(I_array)+.5,{['G= ',num2str(G_tr),' W/m^2','  ',...
        'T= ',num2str(T_tr),' °C']},...
        'HorizontalAlignment','left')
    % legend([num2str(T_tr),'°C'],'Location','Best');
    hold all
    %% relleno curva VI
%     x0r=0;% inicio de x
%     x1r=max(V_array);% fin de x
%     xxr=[x0r x0r V_array(V_array>x0r & V_array<x1r) x1r x1r];
%     yyr=[0 max(I_array) I_array(V_array>x0r & V_array<x1r) 0 0];
%     fill(xxr,yyr,'c');
%     hold all
%% grafica de relleno P real
x=linspace(0,tmpp,max(V_array));
f=@(x) 0*x+cmpp;
y=f(x);
x0=0;
x1=tmpp;
xx=[x0 x0 x(x>x0 & x<x1) x1 x1];
yy=[0 f(x0) y(x>x0 & x<x1) f(x1) 0];
lightGrey = [0.85 0.85 0.85];
fill(xx,yy,lightGrey);
hold all
% pot_real=integral(f,x0,x1);
% pot_ideal=integral(f,xi0,xi1);
global FF pot_real pot_ideal eficiencia
format shortg
areaPV = (largoPV*anchoPV);
Voc = max(V_array);
Isc = max(I_array);
pot_real = tmpp*cmpp;
pot_ideal = Voc*Isc;
FF = pot_real/pot_ideal;
eficiencia = ((Voc*Isc*FF)/(areaPV*G))*10;

text(tmpp/10,cmpp/2,{['\itP_{\rmmp}\rm: ',num2str(pot_real),' Wp'];...
    ['\itV_{\rmoc}\rm: ',num2str(max(V_array)),' V'];...
    ['\itV_{\rmmp}\rm: ',num2str(tmpp),' V'];...
    ['\itI_{\rmsc}\rm: ',num2str(round(max(I_array),3,'significant')),' A'];...
    ['\itI_{\rmmp}\rm: ',num2str(round(cmpp,3,'significant')),' A'];...
    ['FF: ',num2str(round(FF,2,'significant'))];...
    ['Eficiencia: ',num2str(round(eficiencia,3,'significant')),' %']},...
    'Color','k','HorizontalAlignment','left')
hold all
% text(tmpp/10,cmpp/2,{['Potencia real: ',num2str(pot_real),' Wp']},...
%     'Color','k','FontSize',11,'HorizontalAlignment','left')
% text(tmpp/10,cmpp/2.5,{['Potencia ideal: ',num2str(pot_ideal),' Wp']},...
%     'Color','k','FontSize',11,'HorizontalAlignment','left')
% text(tmpp/10,cmpp/3.5,{['Factor de llenado: ',...
%     num2str(FF)]},...
%     'Color','k','FontSize',11,'HorizontalAlignment','left')
% text(tmpp/10,cmpp/5,{['Eficiencia: ',...
%     num2str(eficiencia),' %']},...
%     'Color','k','FontSize',11,'HorizontalAlignment','left')
% plot(tmpp,cmpp,'ro',tmpp,0,'ro',0,cmpp,'ro',0,max(I_array),'ro',...
%     max(V_array),0,'ro')
hold all
plot(tmpp,cmpp,'k*',tmpp,0,'k*',0,cmpp,'k*',0,max(I_array),'k*',...
    max(V_array),0,'k*')
hold all
plot(V_array,I_array,'r','LineWidth',1);
hold all

%% cargo datos medidos
% filename = 'medicionesTesis.xlsx';
% sheet = 'Pmax';
% % V_array_Range = ('F3');
% % A_V_array = (V_array');
% % xlswrite(filename,A_V_array,sheet,V_array_Range)
% 
% P_array_Range = ('I3');
% A_P_array = (P_array');
% xlswrite(filename,A_P_array,sheet,P_array_Range)

end
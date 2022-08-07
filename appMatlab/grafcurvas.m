function [] = grafcurvas(G_tr,T_tr,color,Rs,Rp,n1,n2,Vocn,Iscn,Kv,Ki,...
    Impp,Ns,mg,mVIP,mIT,modelo)
%% entrada de valores de Temp e Irradiancia
Tt = T_tr;
T = Tt + 273.15; % temp. de trabajo
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
Ipv = (Iscn+Ki*(T-Tn))*(G/Gn);
cont = 1;
%%
paso = .1;
for V = 0: paso: Vocn+5
%     syms I
    x = Impp; % valor inicial de I
    Tol = 0.001; % tolerancia o error
    count = 0; % contador de iteraciones
    error = 1; %derivada de la ftn
    %% ecuaciones preliminares
        f1 = @(x) Ipv-Io*(exp((V+x*Rs)/Vt)+exp((V+x*Rs)/((p-1)*Vt))+2)-...
            ((V+x*Rs)/Rp)-x;% Ipv(eq3)
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
    if x < 0 % si hay valores negativos de I se ignorarán
        break
    end
end

%% guardar valores de salida V I P en fichero txt
GT = ({[num2str(round(Tt)),'_',num2str(round(G))]});
GT = char(GT);
results_PIV = [P_array;I_array;I_array];
% modelostr = char(modelo);
filetxt = char({['GrafCurves_PIV_',modelo,'_',GT,'.txt']});
fileID = fopen(filetxt,'w');
fprintf(fileID,'%8s %8s %8s\n','P(W)', 'V(V)', 'I(A)');
fprintf(fileID,'%8.2f %8.2f %8.2f\n',results_PIV);
fclose(fileID);

% global vmarca p
p=P_array;
v=V_array;
% i=I_array;

p = find(p==max(P_array));
vmarca = v(p);
% imarca = i(p);


if (mVIP == 1) && (mIT == 1)% curva VI a temperatura constante
    %% grafico
    figure(mg+1)
    plot(V_array,I_array, color);
%     set(gca,'xtick', 0:2:max(V_array));% valores del eje X
    axis([0 max(V_array+3) 0 max(I_array)+1]);
    xlabel('Tensión (V)');
    ylabel('Corriente (A)');
    title({'Curvas características V-I';
        ['Temperatura: ',num2str(T_tr),' °C']});
    text(4,max(I_array)+.2,{[num2str(G_tr),' W/m^2']},...
        'HorizontalAlignment','left')
    hold all

elseif (mVIP == 0) && (mIT == 1)% curva VP a temperatura constante
    %% grafico
    figure(mg+1)
    plot(V_array,P_array, color);
%     set(gca,'xtick', 0:2:max(V_array));% valores del eje X
    axis([0 max(V_array)+5 0 max(P_array)+50]);
    xlabel('Tensión (V)');
    ylabel('Potencia (W)');
    title({'Curvas características V-P';
        ['Temperatura: ',num2str(T_tr),' °C']});
    text(vmarca, max(P_array)+2.5,{[num2str(G_tr),' W/m^2']},...
        'HorizontalAlignment','left')
    hold all
    %
elseif (mVIP == 1) && (mIT == 0)% curva VI a irradiancia constante
    %% grafico
    figure(mg+1)
    global m xi
    if color(1)=='k'
        m = max(I_array)-0.25*max(I_array);
        xi = 0.25*max(V_array);
    else
        m = m - 0.3;
    end
    plot(V_array,I_array, color);
%     set(gca,'xtick', 0:2:max(V_array));% valores del eje X
    axis([0 max(V_array+3) 0 max(I_array)+3]);
    xlabel('Tensión (V)');
    ylabel('Corriente (A)');
    title({'Curvas características V-I';
        ['Irradiancia: ',num2str(G_tr),' W/m^2']});
    text(xi,m,{[num2str(T_tr),' °C']},'Color',color,...
        'HorizontalAlignment','right')
    hold all
    %
elseif (mVIP == 0) && (mIT == 0)
    figure(mg+1)
    global x y
    if color(1)=='k'
        y = 30;
        x = max(V_array)-5;
    else
        y = y - 4;
    end
    plot(V_array,P_array,'color',color);% curva VP a irradiancia constante
%     set(gca,'xtick', 0:2:max(V_array));% valores del eje X
    axis([0 max(V_array+5) 0 max(P_array)+50]);
    xlabel('Tensión (V)');
    ylabel('Potencia (W)');
    title({'Curvas características V-P';
        ['Irradiancia: ',num2str(G_tr),' W/m^2']});
    text(vmarca,max(P_array)-2,{[num2str(T_tr),' °C']},'Color',color,...
        'HorizontalAlignment','right')
    hold all
end
end
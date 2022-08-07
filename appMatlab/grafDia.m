function [] = grafDia (j,G_array,T_array,Rs,Rp,n1,n2,Vocn,...
    Iscn,Kv,Ki,Impp,Ns)
%% entrada de valores de Temp e Irradiancia
% NOCT = 45;
Tt = T_array;
T = Tt + 273.15; % temp. de trabajo
Tn = 298.15; % temp. de referencia 25+273=298
dT = T-Tn;
G = G_array; % irradiancia de trabajo del panel
Gn = 1000; % irradiancia de referencia 1000 W/m2
%% constantes físicas
q = 1.60e-19 ; % carga del Electrón
K = 1.38e-23 ; % constante de Boltzman
%% condiciones iniciales
p = 2.2;
Vt = (Ns*K*T)/q;
Io = (Iscn+Ki*dT)/(exp((Vocn+Kv*dT)/Vt)-1);% Io1 = Io2
%% ecuaciones preliminares
Ipv = (Iscn+Ki*dT)*(G/Gn);
cont = 1;

global PDia_array

if G ~= 0
    %%
    paso = 1;
    for V = 0: paso: Vocn
%         syms I
        x = Impp; % valor inicial de I
        Tol = 0.001; % tolerancia o error
        count = 0; % contador de iteraciones
        error = 1; %derivada de la ftn
        %% ecuaciones preliminares
        f1 = @(x) Ipv-Io*(exp((V+x*Rs)/Vt)+exp((V+x*Rs)/((p-1)*Vt))+2)-...
            ((V+x*Rs)/Rp)-x;% Ipv(eq3)
        f = f1(x);% evaluo f para el valor inicial de I
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
    %% salvo variable P
    PDia_array(j,1) = max(P_array);
    save ('PDia_array', 'PDia_array')
else
    PDia_array(j,1) = 0;
    save ('PDia_array', 'PDia_array')
end
end
function varargout = FVSim_v2(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @FVSim_v2_OpeningFcn, ...
    'gui_OutputFcn',  @FVSim_v2_OutputFcn, ...
    'gui_LayoutFcn',  [], ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before FVSim_v2 is made visible.
function FVSim_v2_OpeningFcn(hObject, eventdata, handles, varargin)
% poner el GUI en el centro de la pantalla
movegui('center')
global mg mVIP mes dia PDia_array sheet marca_stc
mg = 0;
mVIP = 1;
mes=1;
dia=1;
sheet = 'Meteonorm 6.1';
PDia_array = [];
marca_stc=0;

% Choose default command line output for FVSim_v2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = FVSim_v2_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;


% --- Executes on button press in graf_irrad_btn.
function graf_irrad_btn_Callback(hObject, eventdata, handles)
global color G_init T_init Rs Rp n1 n2 num_curvas Ns Iscn Vocn Impp ...
    Vmpp PmaxE Kv Ki mg mVIP mIT modelo
mIT = 1;% defino temperatura constante
% parametros calculados anteriormente
Rp = str2double(get(handles.Rp_display,'String'));
Rs = str2double(get(handles.Rs_display,'String'));
n1 = 1;
n2 = 1.2;
mg = mg + 1;
% compruebo que los parametros internos estan calculados
if isnan(Rs)
    opc=questdlg({'Antes de graficar debe calcular los parámetros internos',...
        'Quiere calcular los parámetros internos ahora?'},...
        'Error en gráficas','Si','No','No');
    if strcmp(opc,'Si')
        % llamo a la funcion de calcular parametros internos
        calc_btn_Callback(hObject, eventdata, handles);
    end
else
    h = waitbar(0,'Construyendo gráficas, por favor espere...');
    
    color = 'kcrgbmy';
    T_init = str2double(get(handles.temp_constante,'String'));
    G_init = str2double(get(handles.irrad_init,'String'));
    G_incr = str2double(get(handles.irrad_plus,'String'));
    num_curvas = str2double(get(handles.num_curvas_irrad,'String'));
    
    % datos de catalogo
    Ns = str2double(get(handles.Ns_get,'String'));
    Impp = str2double(get(handles.Impp_get,'String'));
    Vmpp = str2double(get(handles.Vmpp_get,'String'));
    PmaxE = str2double(get(handles.pmp_get,'String'));
    Kv = str2double(get(handles.Kv_get,'String'));
    Ki = str2double(get(handles.Ki_get,'String'));
    modelo = (get(handles.modelo_get,'String'));
    
    % arreglo
    ps = str2double(get(handles.cant_serie_FV,'String'));
    pp = str2double(get(handles.cant_paralelo_FV,'String'));
    Iscn = str2double(get(handles.Isc_get,'String'))*pp;
    Vocn = str2double(get(handles.Voc_get,'String'))*ps;
    
    
    for i = 1:num_curvas
        waitbar(i/num_curvas) %waitbar
        grafcurvas (G_init,T_init,color(i),Rs,Rp,n1,n2,Vocn,Iscn,Kv,Ki,...
            Impp,Ns,mg,mVIP,mIT,modelo)
        G_init = G_init + G_incr;
    end
    close(h) % closed the waitbar
end

function irrad_init_Callback(hObject, eventdata, handles)
irrad_init=str2double(get(hObject,'String'));
if isnan(irrad_init) || (irrad_init <= 0)
    errordlg('El valor debe ser numérico y mayor que cero.','ERROR')
    set(handles.irrad_init,'String',200);
end
% --- Executes during object creation, after setting all properties.
function irrad_init_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function irrad_plus_Callback(hObject, eventdata, handles)
irrad_plus=str2double(get(hObject,'String'));
if isnan(irrad_plus) || (irrad_plus <= 0)
    errordlg('El valor debe ser numérico y mayor que cero.','ERROR')
    set(handles.irrad_plus,'String',200);
end

% --- Executes during object creation, after setting all properties.
function irrad_plus_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function temp_init_Callback(hObject, eventdata, handles)
temp_init=str2double(get(hObject,'String'));
if isnan(temp_init)||(temp_init < -30)||(temp_init > 100)
    errordlg('El valor debe estar entre -30 y 100 °C.','ERROR')
    set(handles.temp_init,'String',0);
end

% --- Executes during object creation, after setting all properties.
function temp_init_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function temp_plus_Callback(hObject, eventdata, handles)
temp_plus=str2double(get(hObject,'String'));
if isnan(temp_plus)|| (temp_plus <= 0)
    errordlg('El valor debe ser numérico y mayor que cero.','ERROR')
    set(handles.temp_plus,'String',25);
end

% --- Executes during object creation, after setting all properties.
function temp_plus_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function num_curvas_irrad_Callback(hObject, eventdata, handles)
num_curvas_irrad=str2double(get(hObject,'String'));
if isnan(num_curvas_irrad)||(num_curvas_irrad <= 0)||(num_curvas_irrad > 7)
    errordlg('El valor debe estar entre 1 y 7.','ERROR')
    set(handles.num_curvas_irrad,'String',5);
end

% --- Executes during object creation, after setting all properties.
function num_curvas_irrad_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in reset_irrad_btn.
function reset_irrad_btn_Callback(hObject, eventdata, handles)
ini=char(' ');
set(handles.irrad_plus,'String',ini);
set(handles.irrad_init,'String',ini);
set(handles.num_curvas_irrad,'String',ini);
set(handles.temp_constante,'String',ini);
set(handles.cant_serie_FV,'String',ini);
set(handles.cant_paralelo_FV,'String',ini);

% --- Executes on button press in calc_btn.
function calc_btn_Callback(hObject, eventdata, handles)
% variables de catalogo
global Ns Iscn Vocn Impp Vmpp PmaxE Kv Ki Rs_final Rp_final Io2
Ns = str2double(get(handles.Ns_get,'String'));
Iscn = str2double(get(handles.Isc_get,'String'));
Vocn = str2double(get(handles.Voc_get,'String'));
Impp = str2double(get(handles.Impp_get,'String'));
Vmpp = str2double(get(handles.Vmpp_get,'String'));
PmaxE = str2double(get(handles.pmp_get,'String'));
Kv = str2double(get(handles.Kv_get,'String'));
Ki = str2double(get(handles.Ki_get,'String'));
% llamo a funcion de extracion de parametros (dosdiodos)
dosdiodos
% visualizo los parametros calculados
set(handles.Rs_display, 'string',Rs_final)
set(handles.Rp_display, 'string',Rp_final)
set(handles.nuno_display, 'string',1)
set(handles.ndos_display, 'string',1.2)
set(handles.Io_display, 'string',Io2)
set(handles.Ipv_display, 'string',Iscn)

% --- Executes on button press in reset_param_btn.
function reset_param_btn_Callback(hObject, eventdata, handles)
ini=char(' ');
set(handles.Rs_display,'String',ini);
set(handles.Rp_display,'String',ini);
set(handles.nuno_display,'String',ini);
set(handles.ndos_display,'String',ini);
set(handles.Io_display,'String',ini);
set(handles.Ipv_display,'String',ini);

function modelo_get_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function modelo_get_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function pmp_get_Callback(hObject, eventdata, handles)
pmp_get=str2double(get(hObject,'String'));
if isnan(pmp_get)
    errordlg('El valor de Pmp debe ser numérico y mayor que cero.','ERROR')
    set(handles.pmp_get,'String',0);
end

% --- Executes during object creation, after setting all properties.
function pmp_get_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Voc_get_Callback(hObject, eventdata, handles)
Voc_get=str2double(get(hObject,'String'));
if isnan(Voc_get)
    errordlg('El valor de Voc debe ser numérico y mayor que cero.','ERROR')
    set(handles.Voc_get,'String',0);
end

% --- Executes during object creation, after setting all properties.
function Voc_get_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Vmpp_get_Callback(hObject, eventdata, handles)
Vmpp_get=str2double(get(hObject,'String'));
if isnan(Vmpp_get)
    errordlg('El valor de Vmpp debe ser numérico y mayor que cero.','ERROR')
    set(handles.Vmpp_get,'String',0);
end


% --- Executes during object creation, after setting all properties.
function Vmpp_get_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Impp_get_Callback(hObject, eventdata, handles)
Impp_get=str2double(get(hObject,'String'));
if isnan(Impp_get)
    errordlg('El valor de Impp debe ser numérico y mayor que cero.','ERROR')
    set(handles.Impp_get,'String',0);
end

% --- Executes during object creation, after setting all properties.
function Impp_get_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Isc_get_Callback(hObject, eventdata, handles)
Isc_get=str2double(get(hObject,'String'));
if isnan(Isc_get)
    errordlg('El valor de Isc debe ser numérico y mayor que cero.','ERROR')
    set(handles.Isc_get,'String',0);
end

% --- Executes during object creation, after setting all properties.
function Isc_get_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function fab_get_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function fab_get_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Kv_get_Callback(hObject, eventdata, handles)
Kv_get=str2double(get(hObject,'String'));
if isnan(Kv_get)
    errordlg('El valor de Kv debe ser numérico.','ERROR')
    set(handles.Kv_get,'String',0);
end

% --- Executes during object creation, after setting all properties.
function Kv_get_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Ki_get_Callback(hObject, eventdata, handles)
Ki_get=str2double(get(hObject,'String'));
if isnan(Ki_get)
    errordlg('El valor de Ki debe ser numérico.','ERROR')
    set(handles.Ki_get,'String',0);
end

% --- Executes during object creation, after setting all properties.
function Ki_get_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Np_get_Callback(hObject, eventdata, handles)
Np_get=str2double(get(hObject,'String'));
if isnan(Np_get)
    errordlg('El valor de Np debe ser numérico y mayor que cero.','ERROR')
    set(handles.Np_get,'String',0);
end

% --- Executes during object creation, after setting all properties.
function Np_get_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Ns_get_Callback(hObject, eventdata, handles)
Ns_get=str2double(get(hObject,'String'));
if isnan(Ns_get)
    errordlg('El valor de Ns debe ser numérico y mayor que cero.','ERROR')
    set(handles.Ns_get,'String',0);
end

% --- Executes during object creation, after setting all properties.
function Ns_get_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in reset_ctlg_btn.
function reset_ctlg_btn_Callback(hObject, eventdata, handles)
ini=char(' ');
set(handles.fab_get,'String',ini);
set(handles.modelo_get,'String',ini);
set(handles.pmp_get,'String',ini);
set(handles.Voc_get,'String',ini);
set(handles.Vmpp_get,'String',ini);
set(handles.Impp_get,'String',ini);
set(handles.Isc_get,'String',ini);
set(handles.Kv_get,'String',ini);
set(handles.Ki_get,'String',ini);
set(handles.Ns_get,'String',ini);
set(handles.Np_get,'String',ini);

% --------------------------------------------------------------------
function pi_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function ayuda_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function acercade_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function salir_Callback(hObject, eventdata, handles)
opc=questdlg('Desea salir del programa?','SALIR','Si','No','No');
if strcmp(opc,'No')
    return;
end
clear all
clc,close all

%% funcion de extracion de parametros
function dosdiodos
global Ns Iscn Vocn Impp Vmpp PmaxE Kv Ki Rs_final Rp_final
%entrada de valores de Temp e Irradiancia
G = 1000; % irradiancia de trabajo del panel
Gn = 1000; % irradiancia de referencia 1000 W/m2
Tt = 25;
T = Tt + 273.15; % temp. de trabajo
Tn = 298.15; % temp. de referencia 25+273=298
dT = T-Tn;
% constantes físicas
q = 1.60e-19 ; % carga del Electrón
K = 1.38e-23 ; % constante de Boltzman
%condiciones iniciales
%condiciones iniciales
n1 = 1;
n2 = 1.2;
p = 2.2;
Vt = (Ns*K*T)/q;
Io = (Iscn+Ki*dT)/(exp((Vocn+Kv*dT)/Vt)-1);% Io1 = Io2
Ipv = (Iscn+Ki*dT)*(G/Gn);
TolP = 0.01;
errorP = 1;
marca = 0;
% global marcaPlot
%marcaPlot = 0;
Rs0 = 0.001;
Rs = Rs0;%
Rp0 = (Vmpp/(Iscn-Impp))-((Vocn-Vmpp)/Impp);
Rp = Rp0;%
paso = .1;
cont2 = 1;
errorP_ant=100;

h = waitbar(0,'Calculando parámetros, puede demorar varios minutos.');

while (errorP >= TolP)
    Rs_final = Rs;
    Rp_final = Rp;
    cont = 1;
    for V = 0: paso: Vocn
%         syms I
        x = Impp; % valor inicial de I
        Tol = 0.001; % tolerancia o error
        count = 0; % contador de iteraciones
        error = 1; %derivada de la ftn
        %ecuaciones preliminares
        f1 = @(x) Ipv-Io*(exp((V+x*Rs)/Vt)+exp((V+x*Rs)/((p-1)*Vt))+2)-...
            ((V+x*Rs)/Rp)-x;% Ipv(eq3)
        f = f1(x);% evaluo f para el valor inicial de I
        fini_deriv = diff(sym(f1),1);
        
        fprintf('iter       V           I          f(I)        errorI       P         errorP        Rs          Rp \n')
        while  ((error > Tol || abs(f) > Tol) && x > 0)%
            count = count + 1;
            fprime = eval(fini_deriv);
            xnew = x - (f/fprime);%
            error = abs(x-xnew);%
            x = xnew;
            f = f1(x);% evaluo el nuevo valor de f(I)
            P = V * x;
            fprintf('%3i%12.4f%12.4f%12.4f%12.4f%12.4f%12.4f%12.3f%12.4f \n',...
                count,V,x,f,error,P,errorP,Rs,Rp)
        end
        I_array(cont) = x;
        V_array(cont) = V;
        P_array(cont) = P;
        cont = cont+1;
    end
    PmaxC = max(P_array);
    errorP = abs(PmaxC - PmaxE);% calculo el error
    
    % pasos para incrementar Rs
    if marca == 0
        Rs = Rs + 0.001;
    elseif (errorP<0.01)||(errorP_ant~=100)||(errorP_ant-errorP)<0%si el error aumenta se detiene el programa
        Rs = Rs + 0.001;
        msgbox('Se alcanzó el menor error posible','Fin de programa')
        break
        errorP_ant = errorP;%guardo el valor anterior del error
    elseif errorP > 4
        Rs = Rs + 0.1;
    elseif errorP >= 0.01
        Rs = Rs + 0.001;
    end
    marca = 1;
    Rp = (Vmpp+Impp*Rs)/(Iscn-Io*(exp((Vmpp+Impp*Rs)/Vt)+exp((Vmpp+Impp*Rs)/((p-1)*Vt))+2)-(PmaxE/Vmpp));
    if Rp <= 0
        errordlg('Programa detenido porque el valor de Rp <= 0','ERROR')
        break
    end
    PmaxC_array(cont2) = PmaxC;
    Rs_array(cont2) = Rs;
    Rp_array(cont2) = Rp;
    errorP_array(cont2) = errorP;
    cont2 = cont2 + 1;
    
    waitbar(errorP / TolP) %waitbar
end

close(h) % cierro la barra de progreso

fprintf(' \n')
fprintf('\nPmaxC = %g, errorP = %g, Rs = %g, Rp = %g \n',...
    max(P_array),errorP,Rs_final,Rp_final)
%% fin de la funcion de extraccion de parametros

% --- Executes when selected object is changed in curvas_uibuttongroup.
function curvas_uibuttongroup_SelectionChangedFcn(hObject, eventdata, handles)
global mVIP

if (hObject == handles.PV_radiobtn)
    mVIP = 0;
elseif (hObject == handles.VI_radiobtn)
    mVIP = 1;
end

% --------------------------------------------------------------------
function exportarpi_Callback(hObject, eventdata, handles)
% hObject    handle to exportarpi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function importarpi_Callback(hObject, eventdata, handles)
[File_importpi,PathName_importpi] = uigetfile...
    ('*.dat;*.txt','Seleccione el archivo con los datos');
if isequal(File_importpi,0)
    return
else
    % copiar a la carpeta temp del usuario
    copyfile([PathName_importpi,File_importpi],[tempdir,File_importpi])
    % cargo el pi txt
    filename = fullfile([tempdir,File_importpi]);
    fileID = fopen(filename);
    datosPanelFV = textscan(fileID,'%s %s');
    fclose(fileID);
    fabricante = datosPanelFV{2}{1};
    modeloFV = datosPanelFV{2}{2};
    Rs = datosPanelFV{2}{16};
    Rp = datosPanelFV{2}{17};
    n1 = datosPanelFV{2}{18};
    n2 = datosPanelFV{2}{19};
    Io = datosPanelFV{2}{20};
    Ipv = datosPanelFV{2}{21};
    %
    PmaxE = datosPanelFV{2}{3};
    Voc = datosPanelFV{2}{4};
    Vmpp = datosPanelFV{2}{5};
    Impp = datosPanelFV{2}{6};
    Isc = datosPanelFV{2}{7};
    Ki = datosPanelFV{2}{8};
    Kv = datosPanelFV{2}{9};
    Ns = datosPanelFV{2}{10};
    Np = datosPanelFV{2}{11};
    largoPV = datosPanelFV{2}{12};
    anchoPV = datosPanelFV{2}{13};
    
    set(handles.fab_get,'String',fabricante);
    set(handles.modelo_get,'String',modeloFV);
    
    set(handles.Rs_display,'String',Rs);
    set(handles.Rp_display,'String',Rp);
    set(handles.nuno_display,'String',n1);
    set(handles.ndos_display,'String',n2);
    set(handles.Io_display,'String',Io);
    set(handles.Ipv_display,'String',Ipv);
    %
    set(handles.pmp_get,'String',PmaxE);
    set(handles.Voc_get,'String',Voc);
    set(handles.Vmpp_get,'String',Vmpp);
    set(handles.Impp_get,'String',Impp);
    set(handles.Isc_get,'String',Isc);
    set(handles.Ki_get,'String',Ki);
    set(handles.Kv_get,'String',Kv);
    set(handles.Ns_get,'String',Ns);
    set(handles.Np_get,'String',Np);
    set(handles.largoPV_edit,'String',largoPV);
    set(handles.anchoPV_edit,'String',anchoPV);
end

% --- Executes on button press in graf_temp_btn.
function graf_temp_btn_Callback(hObject, eventdata, handles)
global color G_const T_init Rs Rp n1 n2 num_curvas Ns Iscn Vocn ...
    Impp Vmpp PmaxE Kv Ki T_incr numFig mVIP mIT modelo
mIT = 0;% defino irradiancia constante
% parametros calculados anteriormente
Rp = str2double(get(handles.Rp_display,'String'));
Rs = str2double(get(handles.Rs_display,'String'));
n1 = 1;
n2 = 1.2;
numFig = numFig + 1;
if isnan(Rs)
    opc=questdlg({'Antes de graficar debe calcular los parámetros internos',...
        'Quiere calcular los parámetros internos ahora?'},...
        'Error en gráficas','Si','No','No');
    if strcmp(opc,'Si')
        % llamo a la funcion de calcular parametros internos
        calc_btn_Callback(hObject, eventdata, handles);
    end
else
    h = waitbar(0,'Construyendo gráficas, por favor espere...');
    
    color = 'kcrgbmy';
    %     color = 'o+*.xsd^v<>ph';
    T_init = str2double(get(handles.temp_init,'String'));
    T_incr = str2double(get(handles.temp_plus,'String'));
    G_const = str2double(get(handles.irrad_constante,'String'));
    num_curvas = str2double(get(handles.num_curvas_temp,'String'));
    
    % datos de catalogo
    Ns = str2double(get(handles.Ns_get,'String'));
    Impp = str2double(get(handles.Impp_get,'String'));
    Vmpp = str2double(get(handles.Vmpp_get,'String'));
    PmaxE = str2double(get(handles.pmp_get,'String'));
    Kv = str2double(get(handles.Kv_get,'String'));
    Ki = str2double(get(handles.Ki_get,'String'));
    % arreglo
    ps = str2double(get(handles.cant_serie_FV,'String'));
    pp = str2double(get(handles.cant_paralelo_FV,'String'));
    Iscn = str2double(get(handles.Isc_get,'String'))*pp;
    Vocn = str2double(get(handles.Voc_get,'String'))*ps;
    
    T_init =((T_init+T_incr*num_curvas)-T_incr);% la temp inicia en el mayor valor
    
    for i = 1:num_curvas
        waitbar(i/num_curvas) %waitbar
        % llamo a la funcion de graficar las curvas
        grafcurvas (G_const,T_init,color(i),Rs,Rp,n1,n2,Vocn,Iscn,...
            Kv,Ki,Impp,Ns,numFig,mVIP,mIT,modelo)
        T_init = T_init - T_incr; % decremento la temp para mejor grafica
    end
    close(h) % closed the waitbar
end

function num_curvas_temp_Callback(hObject, eventdata, handles)
num_curvas_temp=str2double(get(hObject,'String'));
if isnan(num_curvas_temp)||(num_curvas_temp <= 0)||(num_curvas_temp > 7)
    errordlg('Este valor debe estar entre 1 y 7.','ERROR')
    set(handles.num_curvas_temp,'String',5);
end

% --- Executes during object creation, after setting all properties.
function num_curvas_temp_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% function reset_temp_btn
function reset_temp_btn_Callback(hObject, eventdata, handles)
ini=char(' ');
set(handles.num_curvas_temp,'String',ini);
set(handles.temp_init,'String',ini);
set(handles.temp_plus,'String',ini);
set(handles.irrad_constante,'String',ini);
set(handles.cant_serie_FV,'String',ini);
set(handles.cant_paralelo_FV,'String',ini);
%%
function irrad_constante_Callback(hObject, eventdata, handles)
irrad_constante=str2double(get(hObject,'String'));
if isnan(irrad_constante) || (irrad_constante <= 0)
    errordlg('Este valor debe ser numérico y mayor que cero.','ERROR')
    set(handles.irrad_constante,'String',1000);
end

% --- Executes during object creation, after setting all properties.
function irrad_constante_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function temp_constante_Callback(hObject, eventdata, handles)
temp_constante=str2double(get(hObject,'String'));
if isnan(temp_constante) || (temp_constante <= -30)
    errordlg('Este valor debe ser numérico.','ERROR')
    set(handles.temp_constante,'String',25);
end

% --- Executes during object creation, after setting all properties.
function temp_constante_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function cant_paralelo_FV_Callback(hObject, eventdata, handles)
cant_paralelo_temp=str2double(get(hObject,'String'));
if isnan(cant_paralelo_temp) || (cant_paralelo_temp <= 0)
        errordlg('Este valor debe ser numérico y mayor que cero.','ERROR')
    set(handles.cant_paralelo_FV,'String',1);
end

% --- Executes during object creation, after setting all properties.
function cant_paralelo_FV_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function cant_serie_FV_Callback(hObject, eventdata, handles)
cant_serie_temp=str2double(get(hObject,'String'));
if isnan(cant_serie_temp) || (cant_serie_temp <= 0)
    errordlg('Este valor debe ser numérico y mayor que cero.','ERROR')
    set(handles.cant_serie_FV,'String',1);
end

% --- Executes during object creation, after setting all properties.
function cant_serie_FV_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in db_popupmenu.
function db_popupmenu_Callback(hObject, eventdata, handles)
global sheet;
sheet = 'Meteonorm 6.1';
val = get(hObject,'Value');
str = get(hObject,'String');
switch str{val}
    case 'Meteonorm 6.1' %
        sheet = 'Meteonorm 6.1';
    case 'NASA-SEE' %
        sheet = 'NASA-SEE';
end

% --- Executes during object creation, after setting all properties.
function db_popupmenu_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in graf_date_btn.
function graf_date_btn_Callback(hObject, eventdata, handles)
global mes dia Rs Rp n1 n2 Vocn Iscn Kv Ki Impp Ns mg
% parametros calculados anteriormente
Rp = str2double(get(handles.Rp_display,'String'));
Rs = str2double(get(handles.Rs_display,'String'));
n1 = 1;
n2 = 1.2;
% datos de catalogo
Ns = str2double(get(handles.Ns_get,'String'));
Impp = str2double(get(handles.Impp_get,'String'));
Kv = str2double(get(handles.Kv_get,'String'));
Ki = str2double(get(handles.Ki_get,'String'));
% arreglo
ps = str2double(get(handles.cant_serie_FV,'String'));
pp = str2double(get(handles.cant_paralelo_FV,'String'));
Iscn = str2double(get(handles.Isc_get,'String'))*pp;
Vocn = str2double(get(handles.Voc_get,'String'))*ps;

mg = mg + 1;
% compruebo que los parametros internos estan calculados
if isnan(Rs)
    opc=questdlg({'Antes de graficar debe calcular los parámetros internos',...
        'Quiere calcular los parámetros internos ahora?'},...
        'Error en gráficas','Si','No','No');
    if strcmp(opc,'Si')
        % llamo a la funcion de calcular parametros internos
        calc_btn_Callback(hObject, eventdata, handles);
    end
else
%
% filename = 'Radiacion_año_moa.xlsx';
% hoja = sheet;
% [num,txt,raw] = xlsread(filename,hoja);% importo del excel
% fecha = raw(2:end,1);% variables
% G_array = [];
% T_array = [];
% fecha_array = [];

fileNASA = fullfile('Clima_Moa_NASA.txt');
fileID = fopen(fileNASA);
datosPanelFV = textscan(fileID,'%s %s %f32 %f32');
fclose(fileID);
% variables
fecha = datosPanelFV{1};
% hora = datosPanelFV{2};
G_array = datosPanelFV{3};
NOCT = 45;
T_array = datosPanelFV{4}+(NOCT-20)*G_array/800;
% fecha_array = [];
% hora_array = 0:1:23;

% comparo las fechas
if dia >= 10 && mes >= 10
    fechaSelect = ([num2str(dia),'/',num2str(mes),'/1900',]);
end

if dia <= 9 && mes <= 9
    fechaSelect = (['0',num2str(dia),'/0',num2str(mes),'/1990',]);
end
if dia <= 9 && mes > 9
    fechaSelect = (['0',num2str(dia),'/',num2str(mes),'/1990',]);
end
if dia > 9 && mes <= 9
    fechaSelect = ([num2str(dia),'/0',num2str(mes),'/1990',]);
end
% comienzo a comparar
for i=1:1:length(fecha)
    iguales = strcmp(fechaSelect,fecha(i,1));% comparo las fechas
    if iguales
        G_array = G_array(i:i+23,1);
        T_array = T_array(i:i+23,1);
%         fecha_array = fecha(i:i+23,1);
        % muestro la tabla
        datostabla = [G_array T_array];
        set(handles.uitable, 'Data', datostabla);
        break
    end
end

h = waitbar(0,'Construyendo la gráfica, puede demorar unos minutos...');
for j=1:length(G_array)
    waitbar(j/24) %waitbar
    grafDia (j,G_array(j),T_array(j),Rs,Rp,n1,n2,Vocn,Iscn,Kv,Ki,Impp,Ns)
end
close(h) % closed the waitbar
% valido la fecha que sea correcta
if ~iguales % si la fecha no aprece muestro este mensaje
     errordlg('La fecha es incorrecta','Error en fecha')
     return
end
% lleno la tabla
load ('PDia_array', 'PDia_array')
datostabla(:,3) = PDia_array;
set(handles.uitable, 'Data', datostabla);
% grafico


figure(mg+1)

% filtar radiacioanes de noche 
for i=1:length(G_array)
    if G_dia(i) >= 1
        G(1,length(G)+1) = G_dia(i);
        time(1,length(G)) = time_dia(i);
    end
end
%fin

%%
figure(mg)
subplot(3,1,1)
% potenciamedidadia22 = dlmread('22marzoEdif21.dat');%borrar
hora_array = linspace(0,1,24);
plot(hora_array,PDia_array/1000,'r','LineWidth',1.2)
xlim([0 23])
dateFormat = 15;
datetick('x',dateFormat)
% hold all
plot(hora_array,PDia_array/1000,'k*','LineWidth',0.5)
hold all%borrar
% plot(hora_array,potenciamedidadia22,'b','LineWidth',0.5)%borrar
xlabel('Horas del dia','FontSize',09);
ylabel('Potencia (W)','FontSize',09);
title({['Producción del panel. Día: ',fechaSelect]},'FontSize',09);
hold all
%
figure(mg)
subplot(3,1,2)
plot(hora_array,G_array,'g')% hora vs G
xlim([0 23])
dateFormat = 15;
datetick('x',dateFormat)
hold all
plot(hora_array,G_array,'k*','LineWidth',0.5)% hora vs G
xlabel('Horas del dia','FontSize',09);
ylabel('Irradiancia (W/m^2)','FontSize',09);
% title({['Día: ',fechaSelect]});
hold all

%
figure(mg)
subplot(3,1,3)
plot(hora_array,T_array,'b')% hora vs T
xlim([0 23])
dateFormat = 15;
datetick('x',dateFormat)
hold all
plot(hora_array,T_array,'k*','LineWidth',0.5)% hora vs T
xlabel('Horas del dia','FontSize',09);
ylabel('Temperatura (°C)','FontSize',09);
% title({['Día: ',fechaSelect]});
end
%%
% --- Executes on selection change in dia_popupmenu.
function dia_popupmenu_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');
str = get(hObject,'String');
global dia
switch str{val}
    case '01'; dia=1;case '02'; dia=2;case '03'; dia=3;case '04'; dia=4;
    case '05'; dia=5;case '06'; dia=6;case '07'; dia=7;case '08'; dia=8;
    case '09'; dia=8;case '10'; dia=10;case '11'; dia=11;case '12'; dia=12;
    case '13'; dia=13;case '14'; dia=14;case '15'; dia=15;case '16'; dia=16;
    case '17'; dia=17;case '18'; dia=18;case '19'; dia=19;case '20'; dia=20;
    case '21'; dia=21;case '22'; dia=22;case '23'; dia=23;case '24'; dia=24;
    case '25'; dia=25;case '26'; dia=26;case '27'; dia=27;case '28'; dia=28;
    case '29'; dia=29;case '30'; dia=30;case '31'; dia=31;
end

% --- Executes during object creation, after setting all properties.
function dia_popupmenu_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in mes_popupmenu.
function mes_popupmenu_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');
str = get(hObject,'String');
global mes
switch str{val}
    case 'enero';mes=1;case 'febrero';mes=2;case 'marzo';mes=3;
    case 'abril';mes=4;case 'mayo';mes=5;case 'junio'; mes=6;
    case 'julio';mes=7;case 'agosto';mes=8;case 'septiembre';mes=9;
    case 'octubre';mes=10;case 'noviembre';mes=11;case 'diciembre';mes=12;
end

% --- Executes during object creation, after setting all properties.
function mes_popupmenu_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in grafRen_btn.
function grafRen_btn_Callback(hObject, eventdata, handles)
global Rs Rp n1 n2 Ns Np Iscn Vocn Impp Vmpp PmaxE Kv Ki mg largoPV ...
    anchoPV FF pot_real pot_ideal eficiencia G_Rend T_Rend marca_stc modelo
% mIT = 1;% defino temperatura constante
% parametros calculados anteriormente
Rp = str2double(get(handles.Rp_display,'String'));
Rs = str2double(get(handles.Rs_display,'String'));
% irrad_init = str2double(get(handles.irrad_init,'String'));
largoPV = str2double(get(handles.largoPV_edit,'String'));
anchoPV = str2double(get(handles.anchoPV_edit,'String'));
n1 = 1;
n2 = 1.2;
mg = mg + 1;
% compruebo que los parametros internos estan calculados
if isnan(Rp)
    opc=questdlg({'Antes de graficar debe calcular los parámetros internos',...
        'Quiere calcular los parámetros internos ahora?'},...
        'Error en gráficas','Si','No','No');
    if strcmp(opc,'Si')
        % llamo a la funcion de calcular parametros internos
        calc_btn_Callback(hObject, eventdata, handles);
    end
      
else

    if marca_stc==1
        T_Rend = 25;
        G_Rend = 1000;
    end
    if marca_stc==0
        G_Rend = str2double(get(handles.irrad_rend,'String'));
        NOCT = 45;
        T_Rend = str2double(get(handles.temp_rend,'String'));
        T_Rend = T_Rend+(NOCT-20)*G_Rend/800;
    end
    
    % datos de catalogo
    Ns = str2double(get(handles.Ns_get,'String'));
    Np = str2double(get(handles.Np_get,'String'));
    Impp = str2double(get(handles.Impp_get,'String'));
    Vmpp = str2double(get(handles.Vmpp_get,'String'));
    PmaxE = str2double(get(handles.pmp_get,'String'));
    Kv = str2double(get(handles.Kv_get,'String'));
    Ki = str2double(get(handles.Ki_get,'String'));
    Iscn = str2double(get(handles.Isc_get,'String'));
    Vocn = str2double(get(handles.Voc_get,'String'));
    modelo = (get(handles.modelo_get,'String'));
    
    grafRend (G_Rend,T_Rend,Rs,Rp,Vocn,Iscn,Kv,Ki,...
        Impp,Ns,Np,mg,largoPV,anchoPV,modelo)
    
    set(handles.pot_ideal_text,'String',pot_ideal);
    set(handles.pot_real_text,'String',pot_real);
    set(handles.FF_text,'String',FF);
    set(handles.eficiencia_text,'String',eficiencia);
    
    marca_stc=0; 
end

function largoPV_edit_Callback(hObject, eventdata, handles)
largoPV_edit = str2double(get(handles.largoPV_edit,'String'));
if isnan(largoPV_edit) || (largoPV_edit <= 0)
    errordlg('El valor debe ser numérico y mayor que cero.','ERROR')
    set(handles.largoPV_edit,'String',1);
end

% --- Executes during object creation, after setting all properties.
function largoPV_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function anchoPV_edit_Callback(hObject, eventdata, handles)
anchoPV_edit = str2double(get(handles.anchoPV_edit,'String'));
if isnan(anchoPV_edit) || (anchoPV_edit <= 0)
    errordlg('El valor debe ser numérico y mayor que cero.','ERROR')
    set(handles.anchoPV_edit,'String',1);
end

% --- Executes during object creation, after setting all properties.
function anchoPV_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function irrad_rend_Callback(hObject, eventdata, handles)
irrad_rend = str2double(get(handles.irrad_rend,'String'));
if isnan(irrad_rend) || (irrad_rend <= 0)
    errordlg('El valor debe ser numérico y mayor que cero.','ERROR')
    set(handles.irrad_rend,'String',1000);
end

% --- Executes during object creation, after setting all properties.
function irrad_rend_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function temp_rend_Callback(hObject, eventdata, handles)
temp_rend = str2double(get(handles.temp_rend,'String'));
if isnan(temp_rend) || (temp_rend <= 0)
    errordlg('Este valor debe ser numérico y mayor que cero.','ERROR')
    set(handles.temp_rend,'String',25);
end
% --- Executes during object creation, after setting all properties.
function temp_rend_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in stc_btn.
function stc_btn_Callback(hObject, eventdata, handles)
global marca_stc
marca_stc=1;
grafRen_btn_Callback(hObject, eventdata, handles)


% --- Executes on button press in ayuda_extract_parameter.
function ayuda_extract_parameter_Callback(hObject, eventdata, handles)
web('https://doi.org/10.1016/j.solener.2011.06.008', '-browser')


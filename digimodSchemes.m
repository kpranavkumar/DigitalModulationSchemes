function varargout = digimodSchemes(varargin)
% DIGIMODSCHEMES MATLAB code for digimodSchemes.fig
%      DIGIMODSCHEMES, by itself, creates a new DIGIMODSCHEMES or raises the existing
%      singleton*.
%
%      H = DIGIMODSCHEMES returns the handle to a new DIGIMODSCHEMES or the handle to
%      the existing singleton*.
%
%      DIGIMODSCHEMES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DIGIMODSCHEMES.M with the given input arguments.
%
%      DIGIMODSCHEMES('Property','Value',...) creates a new DIGIMODSCHEMES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before digimodSchemes_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to digimodSchemes_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help digimodSchemes

% Last Modified by GUIDE v2.5 09-May-2013 13:25:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @digimodSchemes_OpeningFcn, ...
                   'gui_OutputFcn',  @digimodSchemes_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
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


% --- Executes just before digimodSchemes is made visible.
function digimodSchemes_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to digimodSchemes (see VARARGIN)

% Choose default command line output for digimodSchemes
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes digimodSchemes wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = digimodSchemes_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in bpsk.
function bpsk_Callback(hObject, eventdata, handles)
bpskVal = get(handles.bpsk, 'Value');

axes(handles.axes1);
if(get(handles.snrberhold, 'Value'));
    hold all
else
    hold off
end

snr_in = str2num(get(handles.snrval, 'String'));
pe_bpsk = qfunc(sqrt(2*10.^(snr_in/10)));
semilogy(snr_in,pe_bpsk);
xlabel('SNR');
ylabel('BER');

x = ceil(length(snr_in)/2);
y = qfunc(sqrt(2*10.^(x/10)));

text(x,y,...
['\bullet\leftarrow\fontname{times}',...
'bpsk',],'FontSize',14)

data = rand(2000,1)>0.5;    
hmodb = comm.BPSKModulator;
modData1 = step(hmodb, data);
scatterplot(modData1(2:end))


% --- Executes on button press in qpsk.
function qpsk_Callback(hObject, eventdata, handles)
qpskVal = get(handles.qpsk, 'Value');

axes(handles.axes1);
if(get(handles.snrberhold, 'Value'))
    hold all
else
    hold off
end

snr_in = str2num(get(handles.snrval, 'String'));
pe_qpsk = qfunc(sqrt(2*10.^(snr_in/10)));
semilogy(snr_in,pe_qpsk);
xlabel('SNR');
ylabel('BER');

x = ceil(length(snr_in)/2);
y = qfunc(sqrt(2*10.^(x/10)));

text(x,y,...
['\bullet\leftarrow\fontname{times}',...
'qpsk',],'FontSize',14)

data = rand(2000,1)>0.5;
hmod1 = comm.PSKModulator(4,'BitInput',true);
modData1 = step(hmod1, data);
scatterplot(modData1(2:end))

% --- Executes on button press in mpsk.
function mpsk_Callback(hObject, eventdata, handles)
mpskVal = get(handles.mpsk, 'Value');
% if(mpskVal)
%     set(handles.mvalpsk,'Visible','on');
%     set(handles.text6,'Visible','on');
% else
%     set(handles.mvalpsk,'Visible','off');
%     set(handles.text6,'Visible', 'off');
% end
mvalpsk_in = str2num(get(handles.mvalpsk, 'String'));

axes(handles.axes1);
if(get(handles.snrberhold,'Value'))
    hold all
else
    hold off
end

snr_in = str2num(get(handles.snrval,'String'));
pe_mpsk = 2 * qfunc(sqrt(2*10.^(snr_in/10))*sin(pi/mvalpsk_in)) * (1/log2(mvalpsk_in));
semilogy(snr_in,pe_mpsk);
xlabel('SNR');
ylabel('BER');
% title('M-PSK');
format long;
x = ceil(length(snr_in)/2);
y = 2 * qfunc(sqrt(2*10.^(x/10))*sin(pi/mvalpsk_in)) * (1/log2(mvalpsk_in));
text(x,y,...
['\bullet\leftarrow\fontname{times}',...
'm = ',num2str(mvalpsk_in),'-psk'],'FontSize',14)

data = rand(2000*log2(mvalpsk_in),1)>0.5;
hmod1 = comm.PSKModulator(mvalpsk_in,'BitInput',true);
modData1 = step(hmod1, data);
scatterplot(modData1(2:end))


% --- Executes on button press in dpsk.
function dpsk_Callback(hObject, eventdata, handles)
dpskVal = get(handles.dpsk, 'Value');

axes(handles.axes1);
if(get(handles.snrberhold,'Value'))
    hold all
else
    hold off
end


snr_in = str2num(get(handles.snrval,'String'));
pe_dpsk = 0.5 * exp(-10.^(snr_in/10));
semilogy(snr_in,pe_dpsk);
xlabel('SNR');
ylabel('BER');

x = ceil(length(snr_in)/2);
y = 0.5 * exp(-10.^(x/10));

text(x,y,...
['\bullet\leftarrow\fontname{times}',...
'dpsk',],'FontSize',14)

data = rand(2000,1)>0.5;
hmod4 = comm.DPSKModulator(16,'BitInput',true);
modata4= step(hmod4,data);
scatterplot(modata4(2:end))


% --- Executes on button press in evenk.
function evenk_Callback(hObject, eventdata, handles)
evenkVal = get(handles.evenk, 'Value');
% if(evenkVal)
%     set(handles.mvalevenqam,'Visible','on');
%     set(handles.text7,'Visible','on');
% else
%     set(handles.mvalevenqam,'Visible','off');
%     set(handles.text7,'Visible', 'off');
% end

axes(handles.axes1);
if(get(handles.snrberhold,'Value'))
    hold all
else
    hold off
end

mvalevenqam_in = str2num(get(handles.mvalevenqam,'String'));
v = 2^mvalevenqam_in;
snr_in = str2num(get(handles.snrval,'String'));
pe_evenk = 4/mvalevenqam_in * (1 - 1/sqrt(v)) * qfunc(sqrt(3 * mvalevenqam_in * 10.^(snr_in/10) / (v-1)));
semilogy(snr_in,pe_evenk);
xlabel('SNR');
ylabel('BER');

x = ceil(length(snr_in)/2);
y = 4/mvalevenqam_in * (1 - 1/sqrt(v)) * qfunc(sqrt(3 * mvalevenqam_in * 10.^(x/10) / (v-1)));

text(x,y,...
['\bullet\leftarrow\fontname{times}',...
int2str(2^mvalevenqam_in),'-qam',],'FontSize',14)

data = randi([0 3],2000,1);
hmod6= comm.GeneralQAMModulator;
hmod6.Constellation = [1 1i -1 -1i];
modata6= step(hmod6,data);
scatterplot(modata6(2:end))


% --- Executes on button press in oddk.
function oddk_Callback(hObject, eventdata, handles)
oddkVal = get(handles.oddk, 'Value');
% if(oddkVal)
%     set(handles.mvaloddqam,'Visible','on');
%     set(handles.text8,'Visible','on');
% else
%     set(handles.mvaloddqam,'Visible','off');
%     set(handles.text8,'Visible', 'off');
% end

axes(handles.axes1);
if(get(handles.snrberhold,'Value'))
    hold all
else
    hold off
end

mvaloddqam_in = str2num(get(handles.mvaloddqam,'String'));
v = 2^mvaloddqam_in;
snr_in = str2num(get(handles.snrval,'String'));
pe_oddk = 4/mvaloddqam_in * qfunc(sqrt(3 * mvaloddqam_in * 10.^(snr_in/10) / (v-1)));
semilogy(snr_in,pe_oddk);
xlabel('SNR');
ylabel('BER');

x = ceil(length(snr_in)/2);
y = 4/mvaloddqam_in * qfunc(sqrt(3 * mvaloddqam_in * 10.^(x/10) / (v-1)));

text(x,y,...
['\bullet\leftarrow\fontname{times}',...
int2str(2^mvaloddqam_in),'-qam',],'FontSize',14)


% --- Executes on button press in nrcqam.
function nrcqam_Callback(hObject, eventdata, handles)
nrcqamVal = get(handles.nrcqam, 'Value');
% if(nrcqamVal)
%     set(handles.mvalnrqam,'Visible','on');
%     set(handles.text9,'Visible','on');
% else
%     set(handles.mvalnrqam,'Visible','off');
%     set(handles.text9,'Visible', 'off');
% end

axes(handles.axes1);
if(get(handles.snrberhold,'Value'))
    hold all
else
    hold off
end

mvalnrqam_in = str2num(get(handles.mvalnrqam,'String'));
snr_in = str2num(get(handles.snrval,'String'));
v = 2^mvalnrqam_in;
pe_nrqam = (v-1)*qfunc(sqrt(10.^(snr_in/10)));
semilogy(snr_in,pe_nrqam);
xlabel('SNR');
ylabel('BER');

x = ceil(length(snr_in)/2);
y = (v-1)*qfunc(sqrt(10.^(x/10)));

text(x,y,...
['\bullet\leftarrow\fontname{times}',...
int2str(2^mvalnrqam_in),'-qam NR',],'FontSize',14)

data = randi([0 1],96,1);
hModulator = comm.RectangularQAMModulator(8,'BitInput',true);
hModulator.PhaseOffset = pi/4;
modData = step(hModulator, data); 
scatterplot(modData)


% --- Executes on button press in coh.
function coh_Callback(hObject, eventdata, handles)
cohVal = get(handles.coh, 'Value');
% if(cohVal)
%     set(handles.mvalfsk,'Visible','on');
%     set(handles.text23,'Visible','on');
% else
%     set(handles.mvalfsk,'Visible','off');
%     set(handles.text23,'Visible', 'off');
% end

axes(handles.axes1);
if(get(handles.snrberhold,'Value'))
    hold all
else
    hold off
end

mvalfsk_in = str2num(get(handles.mvalfsk,'String'));
snr_in = str2num(get(handles.snrval,'String'));
ps_coh = (mvalfsk_in - 1) * qfunc(sqrt(10.^(snr_in/10)));
semilogy(snr_in,ps_coh);
xlabel('SNR');
ylabel('BER');

x = ceil(length(snr_in)/2);
y = (mvalfsk_in - 1) * qfunc(sqrt(10.^(x/10)));

text(x,y,...
['\bullet\leftarrow\fontname{times}',...
int2str(2^mvalfsk_in),'-fsk coherent',],'FontSize',14)

hmod5= comm.FSKModulator;
data = rand(2000,1)>0.5;
modata5= step(hmod5,data);
scatterplot(modata5(2:end))


% --- Executes on button press in ncoh.
function ncoh_Callback(hObject, eventdata, handles)
ncohVal = get(handles.ncoh, 'Value');
% if(ncohVal)
%     set(handles.mvalfsk,'Visible','on');
%     set(handles.text23,'Visible','on');
% else
%     set(handles.mvalfsk,'Visible','off');
%     set(handles.text23,'Visible', 'off');
% end

axes(handles.axes1);
if(get(handles.snrberhold,'Value'))
    hold all
else
    hold off
end

mvalfsk_in = str2num(get(handles.mvalfsk,'String'));
snr_in = str2num(get(handles.snrval,'String'));
ps_ncoh = (mvalfsk_in - 1)/2 * exp(-0.5*(10.^(snr_in/10)));
semilogy(snr_in,ps_ncoh);
xlabel('SNR');
ylabel('BER');

x = ceil(length(snr_in)/2);
y = (mvalfsk_in - 1)/2 * exp(-0.5*(10.^(x/10)));

text(x,y,...
['\bullet\leftarrow\fontname{times}',...
int2str(2^mvalfsk_in),'-fsk noncoherent',],'FontSize',14)


% --- Executes on button press in disc.
function disc_Callback(hObject, eventdata, handles)
discVal = get(handles.disc, 'Value');
% if(discVal)
%     set(handles.mvalfsk,'Visible','on');
%     set(handles.text23,'Visible','on');
% else
%     set(handles.mvalfsk,'Visible','off');
%     set(handles.text23,'Visible', 'off');
% end

axes(handles.axes1);
if(get(handles.snrberhold,'Value'))
    hold all
else
    hold off
end

mvalfsk_in = str2num(get(handles.mvalfsk,'String'));
snr_in = str2num(get(handles.snrval,'String'));
ps_disc = (1/log2(mvalfsk_in))*((0.5+0.25*(mvalfsk_in/2 + 1)) * (exp(-2*10.^(snr_in/10))));
semilogy(snr_in,ps_disc);
xlabel('SNR');
ylabel('BER');

x = ceil(length(snr_in)/2);
y = (1/log2(mvalfsk_in))*((0.5+0.25*(mvalfsk_in/2 + 1)) * (exp(-2*10.^(x/10))));

text(x,y,...
['\bullet\leftarrow\fontname{times}',...
int2str(2^mvalfsk_in),'-fsk discriminator',],'FontSize',14)



function mvalpsk_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function mvalpsk_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mvalevenqam_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function mvalevenqam_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mvaloddqam_Callback(hObject, eventdata, handles)

function mvaloddqam_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function mvalnrqam_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function mvalnrqam_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in psdhold.
function psdhold_Callback(hObject, eventdata, handles)


% --- Executes on button press in snrberhold.
function snrberhold_Callback(hObject, eventdata, handles)


function qamfs_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function qamfs_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function qamfc_Callback(hObject, eventdata, handles)
dec=[0 0 0 0 0 0 0 1 0 0 1 0 0 0 1 1 0 1 0 0 0 1 0 1 0 1 1 0 0 1 1 1 1 0 0 0 1 0 0 1 1 0 1 0 1 0 1 1 1 1 0 0 1 1 0 1 1 1 1 0 1 1 1 1]; 
dl=length(dec);

%%%% Serial To Parallel and 2-to-4 Level Converter %%%%
    sp=[];
    o1=[];
    o2=[];
    clear i;
    for i=1:4:64;
        sp=[dec(1,i:i+3)];  % Serial to Parallel 4-bit Register
        I=sp(1,1);          % Separation of I and Q bits
        Id=sp(1,2);
        Q=sp(1,3);
        Qd=sp(1,4);
    
        if I==0             % Assigning Amplitude levels for I-channel
            if Id==0
                o1=[o1 -3]; % if input is 00, output is -3
            else
                o1=[o1 -1]; % if input is 01, output is -1
            end
        
        else
            if Id==0
                o1=[o1 1]; % if input is 10, output is 1
            else
                o1=[o1 3]; % if input is 11, output is 3
            end
        end
    
        if Q==0             % Assigning Amplitude levels for Q-channel
            if Qd==0
                o2=[o2 -3]; % if input is 00, output is -3
            else
                o2=[o2 -1]; % if input is 01, output is -1
            end
        
        else
            if Qd==0
                o2=[o2 1]; % if input is 10, output is 1
            else
                o2=[o2 3]; % if input is 11, output is 3
            end
        end
    
    clear sp, clear I, clear Id, clear Q, clear Qd; 
    end


    
    %%%% Oscillator and Balanced Modulator %%%%
    fc=str2num(get(handles.qamfc,'String'));              % Carrier Frequency
    fs=fc*40;            % Sampling Frequency
    t=1:100;             % Duration of Signal
    s=[];
    clear i;
    for i=1:1:16;        % Modulation (multiplication with carrier signals cos and sin)
        Ac=o1(i);
        As=o2(i);
        s1=Ac*cos(2*pi*(fc/fs)*t);
        s2=As*sin(2*pi*(fc/fs)*t);
        s=[s (s1+s2)];
    end
    Fs=str2num(get(handles.qamfs,'String'));;
%     scatterplot(s(2:end));
h = spectrum.welch;
Hpsd= psd(h,s,'Fs',Fs);
axes(handles.axes3);
plot(Hpsd);


% --- Executes during object creation, after setting all properties.
function qamfc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function bfskfs_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function bfskfs_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bfskfc1_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function bfskfc1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bfskfc2_Callback(hObject, eventdata, handles)
bfskfs_in = str2num(get(handles.bfskfs,'String'));
bfskfc1_in = str2num(get(handles.bfskfc1,'String'));
bfskfc2_in = str2num(get(handles.bfskfc2,'String'));

% The number of bits to send - Frame Length
N = 8;
bit_stream = round(rand(1,N));

f1 = bfskfc1_in; 
f2 = bfskfc2_in;
fs = bfskfs_in;

t = 0: 1/fs : 1;

time = [];

FSK_signal = [];
Digital_signal = [];

for ii = 1: 1: length(bit_stream)
    
    % The FSK Signal
    FSK_signal = [FSK_signal (bit_stream(ii)==0)*sin(2*pi*f1*t)+...
        (bit_stream(ii)==1)*sin(2*pi*f2*t)];
    
    % The Original Digital Signal
    Digital_signal = [Digital_signal (bit_stream(ii)==0)*...
        zeros(1,length(t)) + (bit_stream(ii)==1)*ones(1,length(t))];
    
    time = [time t];
    t =  t + 1;
   
end

Fs=10*fs;
h = spectrum.welch;
Hpsd= psd(h,FSK_signal,'Fs',Fs);
axes(handles.axes3);
plot(Hpsd);


% --- Executes during object creation, after setting all properties.
function bfskfc2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function bpskfs_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function bpskfs_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function bpskfc_Callback(hObject, eventdata, handles)
fc = str2num(get(handles.bpskfc, 'String'));
T = 6/fc;
M=2;
Es=T/2;
% fc=6/T;
N=100;
dT=T/(N-1);
t=0:dT:T;

Fs=str2num(get(handles.bpskfs, 'String'));

u0=sqrt(2*Es/T)*cos(2*pi*fc*t);
% scatterplot(u0);
u1=sqrt(2*Es/T)*cos(2*pi*fc*t+2*pi/M);
% scatterplot(u1);

h = spectrum.welch;
Hpsd= psd(h,u0,'Fs',Fs);
axes(handles.axes3);
plot(Hpsd);

% --- Executes during object creation, after setting all properties.
function bpskfc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function snrval_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function snrval_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function mvalfsk_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function mvalfsk_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in aboutproject.
function aboutproject_Callback(hObject, eventdata, handles)
% hObject    handle to aboutproject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
about;

% --- Executes on button press in aboutauthors.
function aboutauthors_Callback(hObject, eventdata, handles)
% hObject    handle to aboutauthors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
authors;


% --- Executes on button press in reset.
function reset_Callback(hObject, eventdata, handles)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.bfskfs, 'String', '100');
set(handles.bfskfc1, 'String', '3');
set(handles.bfskfc2, 'String', '5');
set(handles.qamfs, 'String', '1000');
set(handles.qamfc, 'String', '500');
set(handles.bpskfs, 'String', '1000');
set(handles.bpskfc, 'String', '500');
set(handles.mvalpsk, 'String', '2');
set(handles.mvalevenqam, 'String', '2');
set(handles.mvaloddqam, 'String', '3');
set(handles.mvalnrqam, 'String', '4');
set(handles.mvalfsk, 'String', '2');
cla(handles.axes1);
cla(handles.axes3);

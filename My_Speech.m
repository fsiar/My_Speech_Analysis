function varargout = My_Speech(varargin)
% MY_SPEECH MATLAB code for My_Speech.fig
%      MY_SPEECH, by itself, creates a new MY_SPEECH or raises the existing
%      singleton*.
%
%      H = MY_SPEECH returns the handle to a new MY_SPEECH or the handle to
%      the existing singleton*.
%
%      MY_SPEECH('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MY_SPEECH.M with the given input arguments.
%
%      MY_SPEECH('Property','Value',...) creates a new MY_SPEECH or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before My_Speech_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to My_Speech_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help My_Speech

% Last Modified by GUIDE v2.5 10-Feb-2021 14:24:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @My_Speech_OpeningFcn, ...
                   'gui_OutputFcn',  @My_Speech_OutputFcn, ...
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


% --- Executes just before My_Speech is made visible.
function My_Speech_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to My_Speech (see VARARGIN)

% Choose default command line output for My_Speech
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes My_Speech wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = My_Speech_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in myfft.
function myfft_Callback(hObject, eventdata, handles)
global fx;
global x;



lb=1 ;
up= 256 ;
a= -j*2*pi/256 ;
inva=-a ;
p =0 ; 
n(1:256) =lb:up ;
for k= lb:up 
    nn= n*a  ;
   xx=  nn*k ; 
  xx =  exp( xx) ;
%  p= p+ x(n) ;xx

xx= xx.*x(n) ;
fo(k)= sum(xx) ;
end


 amp=abs(fo) ;
    plot(amp) ; 
% hObject    handle to myfft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in MyRecored.
function MyRecored_Callback(hObject, eventdata, handles)
global fx;
global x;
global mini;
set(handles.error, 'String', ' ');
if(isempty(mini))
set(handles.error, 'String', 'check record time ');
end
if(mini<0)
set(handles.error, 'String', 'check record time ');
end
% mini=str2num(mini);
fx = 10400 ;
%fx = 8000  ; 
Nbits = 16
% disp('please say:.')


% Setup the recording object
my_recorder = audiorecorder(fx, Nbits, 1);
% Record the audio
record(my_recorder, mini) ;
if(my_recorder.TotalSamples==0)
    set(handles.error, 'String', 'Mic not working ');
else
    % Retrieve the sampled recording
    my_voice = getaudiodata(my_recorder)
    % Play the sampled recording
    x=my_voice
    %x=wavrecord(fx * mini,fx);
    axes(handles.axes1);
    play(my_voice);
end




% hObject    handle to MyRecored (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function min_Callback(hObject, eventdata, handles)
global mini;
mini=get(hObject,'String');
mini=str2num(mini);
if(isempty(mini))
set(hObject,'String','0')
end
% hObject    handle to min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of min as text
%        str2double(get(hObject,'String')) returns contents of min as a double


% --- Executes during object creation, after setting all properties.
function min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in opens.
function opens_Callback(hObject, eventdata, handles)
set(handles.error, 'String', ' ');
 global fx ;
 global x;
 fx = 10400 ;
%  global inbits;
[filename, pathname] = uigetfile({'*.WAV;*.wav','sound(*.wav;*.wav)'});
sound = fullfile(pathname, filename); 
[x,fx]= audioread(sound);

% [x,fx]=wavread(sound);
axes(handles.axes1);
plot(x) ;


% hObject    handle to opens (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Myexit.
function Myexit_Callback(hObject, eventdata, handles)
button = questdlg('Do you want to quit the program?', 'Quit the program','Yes','No','No');
switch button
            case 'Yes',
                close('prac1');
            case 'No',
                quit cancel;
end
% hObject    handle to Myexit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

 
% --- Executes on button press in Myexit.
function Splay_Callback(hObject, eventdata, handles)
global fx;
global x;
%  [x,fx]=wavread(sound);
p = audioplayer(x,fx);
playblocking(p);
 
% hObject    handle to Myexit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in engx.
function engx_Callback(hObject, eventdata, handles)
global fx;
global x;
set(handles.msg, 'String', 'Running ');

le = length(x) 
filt = round(le/1024) 
filt =(filt*2)-2 
u=1024 ;
l=1 ;
e(1:filt)= 0 ;   
i=0 ;

 for f=1:filt ;  
     i=i+1 ;
  po =  0 ;
  
 po = sum(x(l:u) .* x(l:u)) ;
 p(i) = sqrt(po) ;
%  e = sqrt(p)  ;
  u =u+512 ;
 l =l+512 ;

 end
 
 axes(handles.axes1); 
    plot(p)
set(handles.msg, 'String', 'finished ');
    
% hObject    handle to engx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in myzero.
function myzero_Callback(hObject, eventdata, handles)
set(handles.msg, 'String', 'Running ');
global fx;
global x;
le = length(x) ;
filt = round(le/1024) ;
filt =(filt*2)-3 
u=1025 ;
l=2 ;
pi(1:filt)= 0 ;   
i=0 ;

 for f=1:filt+1 ;  
     
      z= 0 ;
  for bo=l:u ; 
   
sig = sign(x(bo-1)) * sign(x(bo)) ;
sig2 = sign(x(bo-1)) * sign(x(bo+1)) ;
if( sig == -1 )  ;
    z = z+1 ;
elseif( sig == 0 && sig2 == -1  ) ;
    z = z+1 ;
end

  end
   p(f) = z ;
  u =u+512 ;
  l =l+512 ;
 end
 axes(handles.axes1); 
    plot(p)
set(handles.msg, 'String', 'finished ');

% hObject    handle to myzero (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in spectx.
function spectx_Callback(hObject, eventdata, handles)
 
global fx;
 global fo;
global x;
global specto ;
global spsize  ;
global  winend
 global  startwin

 tic
 L=length(x)
 N=512;
 halfn = N/2
 window=round((L/N) )-3  %tedad filterha
 r1=1;
 r2=N;
 a= -j*2*pi/N ;
 % fo(1:128,1:f) =0 ;
 xx(1:N)=0 ;
 spe(1:N,1:window)= 0 ;

for win=1:window
  co=0  ;
  lb=r1 ;
  up= r2  ;
  n(1:N) =lb:up ;
for k= lb:up 
    co = co+1 ;
    nn= n*a  ;
    xx=  nn*k ; 
    xx =  exp(xx) ;
    
xx= xx*x(n)' ;
f(co)=sum(xx) ;
end 

    r1=r1+N;
    r2=r2+N;  
   %     fo(:,win)= log(f) ;
     spe(1:N,win)= log(f) ;
%      spe(1:N,win)= log(f) ;
%      specto(lb:up)= log(f) ;
end  
spsize= size(spe) 
% fo=imresize(fo,[128 256]);
axes(handles.axes1);

imshow(spe)
% amp=abs(fo) ;
%     plot(amp) ; 
toc
 
% hObject    handle to spectx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in orgx.
function orgx_Callback(hObject, eventdata, handles)
 global fx ;
 global x;

axes(handles.axes1);
plot(x(1:256)) ;
% hObject    handle to orgx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in myhist.
function myhist_Callback(hObject, eventdata, handles)
global fx ;
 global x;
%   global inbits 
  data=x ;
  Fs = fx ;
z=0 ;
le = length(x) 
[mm,n]=size(data) ;
dmax=max(data)      ;     
dmin=abs(min(data))      
dhigh=max(dmax,dmin)   
dismin=min(data)    
dismin=min(dismin)     
dist=dhigh + dmin ;

dista = dist/10
z1=0 ;
z2=0 ;
z3=0 ;
z4=0 ;
z5=0 ;
z6=0 ;
z7=0 ;
z8=0 ;
z9=0 ;
z10=0 ;

dismin+dista ;
for i=1:le
   
    if (dismin+dista >x(i) & x(i)>dismin )
       
        z1= z1+1  ;
    
    elseif (dismin+dista>x(i)  & x(i)<dismin+(2*dista))
        z2= z2+1  ;
     
    elseif ((2*dista)+dismin<x(i)  & x(i)<dismin+(3*dista))
        z3= z3+1  ;
    
    elseif ((3*dista)+dismin<x(i) & x(i)<dismin+(4*dista))
        z4= z4+1 ;
        elseif ((4*dista)+dismin<x(i) & x(i)<dismin+(5*dista))
        z5= z5+1;
        elseif ((5*dista)+dismin<x(i) & x(i)<dismin+(6*dista))
        z6= z6+1;
        elseif ((6*dista)+dismin<x(i) & x(i)<dismin+(7*dista))
        z7= z7+1;
        elseif ((7*dista)+dismin<x(i) & x(i)<dismin+(8*dista))
        z8= z8+1;
        elseif ((8*dista)+dismin<x(i) & x(i)<dismin+(9*dista))
        z9= z9+1;
        elseif ((9*dista)+dismin<x(i) & x(i)<dismin+(10*dista))
        z10= z10+1 ;
    end
end
  m(1:10) = [z1,z2,z3,z4,z5,z6,z7,z8,z9,z10]  
  a=dhigh(1) + abs(dismin(1))
  aa=a/10
  px=dismin(1):aa:dhigh(1)
plot( px(1:10),m ,   '-')
% 
% j= dismin-1:dista:dhigh+1 ;
% 
%     mmax=max(m)  
% 
% plot( dismin(1),0:10:z1 ,'.') ;
% hold on
% 
%   plot( dismin(1)+dista(1),0:10:z2 ,'.') ;
%   plot([dismin(1)+(2*dista(1)),0:10:z3],'.') ;
% 
% 
%     plot( dismin(1)+(3*dista(1)),0:10:z4 ,'.') ;
%       plot( dismin(1)+(4*dista(1)),0:10:z5 ,'.') ;
%        plot( dismin(1)+(5*dista(1)),0:10:z6 ,'.') ;
%     plot( dismin(1)+(6*dista(1)),0:10:z7 ,'.') ;
%       plot(  dismin(1)+(7*dista(1)),0:10:z8 ,'.') ; 
%         plot( dismin(1)+(8*dista(1)),0:10:z9 ,'.') ;
%     plot( dismin(1)+(9*dista(1)),0:10:z10 ,'.') ;
%         plot( dismin(1)+(10*dista(1)),0:10:z10 ,'.') ;
%   
% 
% %  plot(j,z3) ;
% % hold on
% %  plot(j,z4) ;
%  
% xlabel('amplitude')
% ylabel('number')

   hold off
% hObject    handle to myhist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in autcor.
function autcor_Callback(hObject, eventdata, handles)

global fx
global x


L=length(x);

N=1024;
f=round((L/1024)*2 )-4  %%%%%%%%%%%%%%????? ???????%%%%%%%%%%%%%%%5

r1=1;
r2=1024;

acoro(1:1024,1:f) =0 ;  %%%%%%%%%%???????????%%%%%%%%%%%%5
out2(1:L)=0 ;
for win=1:f
    
lb=r1 ;
up= r2 ; 

 co=0 ;
n(1:1024) =1:1024 ;

for k= 1:1024
      xx(1:1024)=0 ;
  xx(1:1024)=lb+n+k ;
  co=co+1  ; 

  %%%%%%%%%%%%%????????%%%%%%%%%%%%%%%%%%%%%%%%%5
%  sum1= sum(x(n+lb)) ;
% sum2= sum(x(xx)) ;
%  sum12 = sum1 *sum2 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
m= sum(x(n+lb) )*sum( x(xx) );
out(co)= sum(x(n+lb) .* x(xx) )/m;

end 

  acoro(1:1024,win)= out(1:1024) ;
  out2(lb:up) = out;
    r1=r1+512;
    r2=r2+512;      
end

[pich,xne] = max(out2(2:L))
   
xne= xne-1 ;
fpich= pich/xne
st = num2str(fpich) 
set(handles.error, 'String', st);

acoro=imresize(acoro,[128 128]);
imshow(acoro)
% hObject    handle to autcor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in hamwin.
function hamwin_Callback(hObject, eventdata, handles)


global fx;
global x;
global inp ;

%set(handles.error, 'String', 'please inter numbe big than 513');

inp=514
h=linspace(-5,5,1024);
y = abs(sinc(h));
L=length(x);

N=1024;
% f=round((L/1024) )  %tedad filterha

a= -j*2*pi/1024 ;
lb= inp-512 
up=inp+511 
n(1:1024) =lb:up ;
 si(1:1024) =1:1024  ;
 xn(1:1024)=0 ;
 newx= x(n) * y(si)' ;
for k= 1:1024
%           xn(si)=y(si)*x(k+lb); 

    nn= n*a  ;
   xx=  nn*k ; 
  xx =  exp( xx) ;
xx= xx .*newx(n)' ;
fo(k)= sum(xx) ;
end 


plot(abs(fo)) ;

set(handles.erorr, 'String', 'output');
% hObject    handle to haming (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% hObject    handle to hamwin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in cap.
function cap_Callback(hObject, eventdata, handles)

global fx;
global x;


tic   %% time estimate 
L=length(x)

N=1024;
f=round((L/1024) )-3  %tedad filterha
                                                                                                                             
r1=1;
r2=1024;
cont=0 ;
fori(1:L) = 0 ;
a= -j*2*pi/1024 ;
inva= -a ;
fo(1:1024,1:f) =0 ;
xx(1:1024)=0 ;
for win=1:f
    
lb=r1 ;
up= r2  ;
n(1:1024) =lb:up ;

for k= lb:up 
    cont=cont+1 ;
    nn= n*a  ;
   xx=  nn*k ; 
   xx =  exp(xx) ;

xx= xx.*x(n) ;
 xx= sum(xx) ;

fori(cont)= xx ;
end 
fori=log(fori) ;

for indft= lb:up 
    
    nn2= n*inva  ;
   xx2=  nn2*indft ; 
   xx2 =  exp(xx2) ;

xx2= xx2.*fori(n) ;


fo(1:1024,win)= (1/1024)*(sum(xx2)) ;

end 

    r1=r1+512;
    r2=r2+512;      
end  

fo=imresize(fo,[128 256]);
imshow(fo)
% amp=abs(fo) ;
%     plot(amp) ; 

toc

% hObject    handle to cap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

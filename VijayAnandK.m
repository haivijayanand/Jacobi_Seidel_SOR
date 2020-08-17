%Assignment -1
%K.VIJAY ANAND (05775) , ME AERO

clear;
close all;
clc;


% SAMPLE INPUT
%  A=[7 3 -1 2; 3 8 1 -4; -1 1 4 -1; 2 -4 -1 6];
%  b= [-1;0;-3;1];
%  x0=[0;0;0;0];
%  IterMax=1000;
%  Tol=1e-10;


%**************************************************************************

disp('********** SOLUTION TO SET OF LINEAR EQUATIONS **********')
disp(' ');
disp('Choose The Method . . .')
disp(' ');
disp('JACOBI METHOD      . . . 1');
disp('SEIDEL METHOD      . . . 2');
disp('RELAXATION METHOD  . . . 3');

disp(' ');
disp('Press Enter to select the Default Option in bracket');
disp(' ');
Reply= input('Select The Method     [1]   : ', 's');
if isempty(Reply)
    Method = '1';         % Default Value
else
    Method=num2str(Reply);
end
disp(Method);
Method=str2num(Method);

%**************************************************************************

[filename, pathname] = uigetfile('*.dat;*.txt', 'Choose the text file having Augmented matrix');
if isequal(filename,0)||isequal(pathname,0)
    disp('File not found !!! . . . ')
    disp('Loading Default File. . . ');
    Ab=load('Matrix_4.dat')

else
    disp([pathname, filename, ' found . . .'])
    Ab=load([pathname,filename])
end

%**************************************************************************



[N c]=size(Ab);
if N+1~=c
    disp('The input format is not in the Augmented Matrix  form !!!');
    disp(' Loading Default Matrix_4.dat . . . ');
    Ab=load('Matrix_4.dat');
end



A=Ab(:,1:end-1)
b=Ab(:,end)

disp('Press Enter');
pause;

%**************************************************************************
clc;

disp('Checking Matrix for Diagonal Dominance . . .')
% This is based on finding for how many rows the condition
% the absolute value of the diagonal element in a row is
% strictly greater than than the sum of absolute value
% of the rest of the elements in that row.

%size gives how many rows and columns in the A matrix
rowcol=size(A);
n=rowcol(1);
% count = for how many rows is the inequality met that
% the absolute value of the diagonal element in a row is
% strictly greater than than the sum of absolute value
% of the rest of the elements in that row
count=0;
for i=1:1:n
    sumrow=0;
    for j=1:1:n
        if i~=j
            sumrow=sumrow+abs(A(i,j));
        end
    end
    if abs(A(i,i))>sumrow
        count=count+1;
    end
end

if count==n
    disp('Matrix is strictly diagonal dominant')
else
    disp('**************************************** ');
    disp(' ');
    disp('Matrix is NOT strictly diagonal dominant')
    disp('Solution is not guaranteed!!!')
    disp(' ');
    disp('**************************************** ');
end

disp('Press Enter');
pause;

%**************************************************************************

clc;

disp(['Enter ' num2str(N) ' Initial Conditions (x0)']);
disp(' ');


for(i=1:N)
    str=['Enter  ' num2str(i)  '  x0 value   [0]    :'];
    Reply=input(str,'s');
    if isempty(Reply)
        x0(1,i) = 0;         % Default Value
    else
        x0(1,i)=num2str(Reply);
    end
end

disp('Initial Conditions . . .')
x0

disp('Press Enter');
pause;
%**************************************************************************

Reply= input('Maximum Number of Permissible Iterations   [1000]   : ', 's');
if isempty(Reply)
    Max_Iter = 1000;         % Default Value
else
    Max_Iter=str2num(Reply);
end

disp(['Maximum Iterations =  ' num2str(Max_Iter)])

%**************************************************************************

disp(' ');

Reply= input('Maximum Error Tolerance Required   [1e-10]   : ', 's');
if isempty(Reply)
    Tol = 1e-10;            % Default Value
else
    Tol=str2num(Reply);
end
disp(['Maximum Tolerance =  ' num2str(Tol)])
disp('Press Enter');
pause;

%**************************************************************************

disp(' ');
if Method==3
    Reply= input('Relaxation Factor   [1.5]    : ', 's');
    if isempty(Reply)
        w = 1.5;            % Default Value
    else
        w=str2num(Reply);
    end
    disp(['Relaxation Factor =  ' num2str(w)])
end

%**************************************************************************



disp(' ');
disp('Solution  . . .');

if Method==1
    [X,Niter,Xiter,Error]= Jacobi2(A,b,x0,Max_Iter,Tol);
elseif Method==2
    [X,Niter,Xiter,Error]= Seidel2(A,b,x0,Max_Iter,Tol);
elseif Method==3
    [X,Niter,Xiter,Error]= Relax2(A,b,x0,w,Max_Iter,Tol);
end

%**************************************************************************

figure(1);

for(i=1:N)
    subplot(N,1,i);
    if i==1 title('Solution History');
        hold on;
    end

    plot(Xiter(i,:));
    grid on;
    ylabel(['x ' num2str(i)]);
    hold on;

end

xlabel('Iterations --->');


figure(2)
loglog(Error);
grid on;
xlabel('Iterations --->');
ylabel('Norm (Error)');
title ('Error History');

if(Niter>=Max_Iter)
    disp('Maximum Iterations Reached or Solution not Converged!!!');
    disp(' ');
    disp('Possible Reasons . . .');
    disp(' 1. Number of iterations not sufficient to converge on tolerance');
    disp(' ');
    disp(' 2. The Matrix is not diagonally dominant');
    close all;
else
    disp('Solution COnverged Successfully on tolerance!!!');
end


    


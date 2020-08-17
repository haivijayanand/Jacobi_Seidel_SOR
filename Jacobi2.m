function [X,N_Iter,Xiter,Error]= Jacobi2(A,b,x0,IterMax,Tol)
%Assignment -1
%K.VIJAY ANAND (05775) , ME AERO
% Solution to set of linear euations (A x = b) by
% Jacobi Method
% *** input ***
% A - Matrix, 
% b - Matrix, 
% x0 -Initial Condition
% w - Relaxation Parameter
% IterMax - Maximum Number of Iterations
% Tol - Error Tolerance
% 
% *** Output ***
% X - Converged Solution
% N - Number of iterations
% Xiter - Solution History
% Error - Error History


% SAMPLE INPUT
%  A=[7 3 -1 2; 3 8 1 -4; -1 1 4 -1; 2 -4 -1 6];    
%  b= [-1;0;-3;1];
%  x0=[0;0;0;0];
%  IterMax=1000;
%  Tol=1e-10;



N=length(A);

x=zeros(N,IterMax);     %defining storing variable for solution at various iterations
Error=zeros(N,1);
x(:,1)=x0;

disp(' ');
disp(' ');
disp('*** SOLUTION BY JACOBI METHOD ***');

disp('Diagonal Matrix');
D=diag(diag(A))



disp('Off Diagonal Matrix');
Aoff=A-D


disp('Inverse of Diagonal Matrix');
Dinv=inv(D)
disp('Press Enter . . .');
pause;

% loop for iterations

for k=1:IterMax
    
    x(:,k+1)=Dinv*(b-Aoff*x(:,k));               %Computes Solution for next iteration

    Error(k)= norm(x(:,k+1)-x(:,k));              %Computes the Error
    
    if Error(k)<Tol
        break
    end
end




x(:,k+2:end)=[];
Error(k+2:end)=[];


X=x(:,end);
Xiter=x;
N_Iter=k;


% show the final solution
disp('Final Solution');
X
% show the total iteration
disp('Number of Iterations');
N_Iter
end
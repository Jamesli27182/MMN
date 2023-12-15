%% part2 try3
clear all;
close all;
clc;
% Data generation 
N_of_ser = 3
N = 6000+1
lambda = 4;
miu = 1/4;
arrive_time_int = exprnd(1/lambda)*ones(1,N)
for i = 2:N;
 arrive_time_int(i) = arrive_time_int(i-1)+exprnd(1/lambda);
end
ser_time = exprnd(miu,1,N)
ser_wind_depart = zeros(N_of_ser,N)
Infor_matrix = [arrive_time_int;
    zeros(1,N);
    ser_time;
    zeros(1,N);
    ser_wind_depart;]
% Infor_matrix(4:(4+N_of_ser-1),1)
%% start to simulate
n_arr = 1
while n_arr <= N-1;
arrive_time = Infor_matrix(1,n_arr);
window_finished_time = Infor_matrix(4:(4+N_of_ser-1),n_arr);
[if_go_to_window,window_num] = ser_or_wait(arrive_time,window_finished_time);
if if_go_to_window == 0;%
      [num,pos]= find_minnum_position(Infor_matrix(4:(4+N_of_ser-1),n_arr));
% min number + server time
      Infor_matrix((4+pos-1),n_arr)=num+Infor_matrix(3,n_arr);
% serer time shift to next coloumu
    Infor_matrix((4:4+N_of_ser-1),(n_arr+1)) = Infor_matrix((4:4+N_of_ser-1),(n_arr));
    Infor_matrix(end,n_arr) = pos;
else
% server finished time = arrive time + server time
Infor_matrix(window_num+4-1,n_arr) = arrive_time + Infor_matrix(3,n_arr);
% serer time shift to next columu
Infor_matrix(4:(4+N_of_ser-1),(n_arr + 1)) = Infor_matrix(4:(4+N_of_ser-1),n_arr);
Infor_matrix(end,n_arr) = window_num;
end

 n_arr = n_arr + 1
end
time_line = Infor_matrix(1,1:(end-1));% should record arrive tiem and depart time 
for i = 1:(N-1);
    time_line = [time_line,Infor_matrix(4+(Infor_matrix(end,i))-1,i)];
end
for i = 1:(length(time_line)-1);
    for j = 1:(length(time_line)-i-1);
        if time_line(j)>time_line(j+1);
            a = time_line(j);
            b = time_line(j+1);
            time_line(j) = b;
            time_line(j+1) = a;
        end
    end
end
k = 0;
quene_line = [];
for i = 1:length(time_line);
    k = k +  if_arr_depart(time_line(i),Infor_matrix,N);
    quene_line = [quene_line,k]  ;
end
plot(time_line,quene_line,'*-')
grid on
A = 0;
for i = 2:length(quene_line);
A = A + (time_line(i)-time_line(i-1))*quene_line(i-1);
end
aver_quene = A/(time_line(2*(N-1)))
%% function used
function [C,window_num] = ser_or_wait(arrive_time,A);% A is window's line coloum n_ser_window x 1
    L = length(A);   % Number of servers 
    window_num = [];
    for i = 1:L;
        if arrive_time>=A(i)
            C = true; % 1 means go to window
            window_num = i;
            break
        else
            C =false; % 0 means go to quene
        end
    end
end 

function [min_n,position] = find_minnum_position(A)
L = length(A);
min_n = min(A);
position = [];
for i = 1:L
    if min_n == A(i)
        position = i;
        break
    end
end
end

function [F] = if_arr_depart(b,Infor_matrix,N)
F = -1;
for i = 1:N
    if b == Infor_matrix(1,i)
        F = 1;
    end
end
end
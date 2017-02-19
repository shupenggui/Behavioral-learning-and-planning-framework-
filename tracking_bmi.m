function [ obj ] = tracking_bmi( person,range,a,D,shift )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    vbmi = importdata(strcat(num2str(person),'vbmi.mat'));
    Dc = importdata(strcat(num2str(person),'D.mat'));
    A = Dc(:,1:3);
    B = Dc(:,4:end-1);
    C = Dc(:,end);
    obj = zeros(range,1);
    move = [vbmi(3+shift);vbmi(2+shift);vbmi(1+shift)];
    obj(1:3) = move;
    for i = 1:range-3
        move = A*move + B * [D*a(:,i+3);D*a(:,i+2);D*a(:,i+1)] + C;
        obj(i+3) = move(1);
    end

end


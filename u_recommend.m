function [ u,a,D ] = u_recommend( s,person,lambda,th,range ,c,dtarget,shift)

%     v = load(s);
%     U = v(:,3:6);
%     U(U(:,1)==-1,:) = [];
%     U(U(:,2)==-1,:) = [];
%     U(U(:,3)==-1,:) = [];
%     U(U(:,4)==-1,:) = [];
    [d,v,nv] = raw_data_separation(s);

    %v = load('orig139.tsv');

    U = nv(:,3:6);temp = v(:,3:6);
    U(temp(:,1)==-1,:) = [];temp(temp(:,1)==-1,:) = [];
    U(temp(:,2)==-1,:) = [];temp(temp(:,2)==-1,:) = [];
    U(temp(:,3)==-1,:) = [];temp(temp(:,3)==-1,:) = [];
    U(temp(:,4)==-1,:) = [];temp(temp(:,4)==-1,:) = [];
    N = size(U,1);
    U = U';

    Dc = importdata(strcat(num2str(person),'D.mat'));
    ut = importdata(strcat(num2str(person),'u.mat'));
    vbmi = importdata(strcat(num2str(person),'vbmi.mat'));
    
    A = Dc(:,1:3);
    B = Dc(:,4:end-1);
    C = Dc(:,end);

%     N = 10000;
%     U = mapminmax(U(1:N,:)');

    index = importdata('index3372_1_6.mat');
    D = U(:,index(1:th,1));
    target = vbmi(1+shift) - dtarget;

    cvx_begin
        variable a(th,range);
        ieq = [vbmi(3+shift);vbmi(2+shift);vbmi(1+shift)];
        for i = 1:range-3
            ieq = A*ieq + B * [D*a(:,i+3);D*a(:,i+2);D*a(:,i+1)] + C;
        end
        minimize(sum(c'*D*a) + lambda * norm(D*(a(:,2:end)-a(:,1:end-1)),'fro'));
       % minimize(sum(c'*D*a) + lambda * sum(sum(abs(D*(a(:,2:end)-a(:,1:end-1))))));
        subject to
            ieq(1) <= target;
            D*a(:,1) == ut(1+shift,:)';
            a >= 0;
            a <= 1;
            ones(1,th) * a == ones(1,range);
%            D(2,:)*a <= 0.5;
%            D(2,:)*a >= 0;
    cvx_end
    
    u = D*a;


end


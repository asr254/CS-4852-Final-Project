clear

%set number of agents
numAgents = 100; %original does not work at n = 2 or n = 3 @x=1, minimize successes by having 4 agents and cars (increases <4 and >4); 1000/1000 @ 2 agents and cars for all x
v = numAgents;

%set number of tested trials
numTrials = 1000;

perc = [];
avgAmt = [];
for i = 1:numTrials
    %generates random permutation of priority order
PO = randperm(v);

%generates random customer preferences
CP = [];
for i = 1:v
    CP = [CP, randperm(v)];
end

%generates random dealer preferences
DP = [];
for i = 1:v
    DP = [DP, randperm(v)];
end

    %perc = [perc, ANDIFcomp(PO,CP,DP)];
    avgAmt = [avgAmt, ANDIFcomp2(PO,CP,DP)];
    
end

%fprintf('Number of successes with %i customers and cars: %i/%i trials\n',numAgents,numTrials*mean(perc),numTrials);
fprintf('Average amount ANDIF(M2) is greater than ANDIF(M1): %.6f\n',mean(avgAmt));

%takes in consumer preferences C, dealer preferences D, and matching M to
%generate the ANDI of M
function a = ANDIF(C,D,M)
    v = sqrt(length(C));
    array = [];
    for k = 1:v
        c = M(k*2-1);
        d = M(k*2);
        
        cPref = C((((k-1)*v)+1):(k*v));
        dPref = D((((d-1)*v)+1):(d*v));
        
        value1 = find(cPref == d);
        value2 = find(dPref == c);
        
        %x takes value [0,1]; increasing x increases weight of customer satisfaction vs dealer satisfaction ----> determined empirically
        x = .5;
        normVal = ((1/value1)^(1-x))*((1/value2)^x);
        
        array = [array, normVal];
    end
    a = mean(array);
end

%takes in priority order P and consumer preferences C to generate matching
%M1 using serial dictatorship algorithm
function M = M1(P,C)
CP = C;
v = sqrt(length(CP));
M = [];
    for i = 1:v
        j = P(i); %j is consumer in priority order
        m = CP(((j-1)*(v+1-i))+1); %m is j's matched car
        p = [j m]; %p is a pair matching priority consumer with favorite remaining car
        M = [M, p]; %add p to matching M
        CP = CP(CP~=m); %remove all instances of m from CP
    end
end

%takes in consumer preferences C and dealer preferences D to generate matching
%M2 using dealer-proposing deferred acceptance algorithm
function M = M2(C,D)
CP = C;
DP = D;

%initialize data structure T such that T{c,1} = matrix of preferences for 
%person c without people who have denied them, and T{c,2} = who they are 
%matched to, or 0 if unmatched
for j = 1:sqrt(length(CP))
    T{j,1} = CP(((j-1)*sqrt(length(CP))+1):((j*sqrt(length(CP)))));
    T{j,2} = 0;
end

%initialize data structure S such that S{c} = matrix of preferences for
%dealer d
for j = 1:sqrt(length(DP))
    S{j} = DP(((j-1)*sqrt(length(DP))+1):((j*sqrt(length(DP)))));
end

while hasZero(T)~=0 %loop of algorithm
    cNew = hasZero(T);
    dPref = T{hasZero(T),1}; %get free person's dealer preferences
    d = topNotProp(dPref); %find free person's top choice of those not proposed to yet
    cOld = hasValue(d,T); %find person who has already proposed to d, or 0 if d is free
    
    if hasValue(d,T) == 0 %if no one has proposed to d, assign c to d
        T{cNew,2} = d;
    else %if someone has proposed to d, compare d's preferences to see if 
        if dCompare(cOld,cNew,S{d}) == 1 %if d prefers cNew, free cOld, mark cOld as having proposed to d, and assign d to cNew
            T{cOld,2} = 0;
            T = haveProposed(cOld,d,T);
            T{cNew,2} = d;
        else %if d prefers cOld, mark cNew as having proposed to d
            T = haveProposed(cNew,d,T);
        end
    end
end

%create matching M from data structure T
M = [];
for r = 1:sqrt(length(CP))
    M = [M,r];
    M = [M,T{r,2}];
end

end

%mark c as having proposed to d by changing its T{c,1} matrix to have value
%d at index k changed to 0
function t = haveProposed(c,d,T)
    t = T;
    temp = t{c,1};
    k = find(temp == d);
    temp(k) = 0;
    t{c,1} = temp;
end

%returns 1 if d prefers cNew to cOld, and 0 otherwise
function b = dCompare(cOld,cNew,cPref)
    b = 0;
    if find(cPref == cOld) > find(cPref == cNew)
        b = 1;
    end
end

%returns a consumer who is free or 0 if there are no more free people
function c = hasZero(T)
    K = size(T);
    max = K(1);
    c = 0;
    for i = 1:max
        if T{i,2} == 0
            c = i;
        end
    end
end

%takes in dealer number and T matrix and returns the consumer who already
%has that dealer assigned, or returns 0 if no consumer has that dealer
function c = hasValue(d,T)
    K = size(T);
    max = K(1);
    c = 0;
    for i = 1:max
        if T{i,2} == d
            c = i;
            break;
        end
    end
end

%takes in dealer preferences for consumer and returns first nonzero number
function x = topNotProp(dpref)
    for i = 1:length(dpref)
        if dpref(i) ~= 0
            x = dpref(i);
            break;
        end
    end
end

%%view properties of matchings
function viewProp(P,C,D)
        fprintf('Priority order:');
        display(P);
        fprintf('Consumer Preferences:');
        display(C);
        fprintf('Dealer Preferences:');
        display(D);
        fprintf('M1:');
        display(M1(P,C));
        display(ANDIF(C,D,M1(P,C)));
        fprintf('M2:');
        display(M2(C,D));
        display(ANDIF(C,D,M2(C,D)));
end

%takes in priority order P, consumer preferences C, and dealer preferences D to
%print true if the ANDI of M2 is greater than or equal to the ANDI of M1; 
%if the value is false, displays the priority order, consumer preferences, and dealer preferences
function a = ANDIFcomp(P,C,D)
    a = 0;
    if ANDIF(C,D,M2(C,D)) >= ANDIF(C,D,M1(P,C))
        a = 1;
    else
        %viewProp(P,C,D);
    end
end

%same as ANDIFcomp but returns the difference in the ANDIFs
function a = ANDIFcomp2(P,C,D)
    a = ANDIF(C,D,M2(C,D)) - ANDIF(C,D,M1(P,C));
end

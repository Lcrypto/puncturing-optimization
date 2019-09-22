%
%Copyright(c) 2012, German Svistunov
%All rights reserved.
%Redistribution and use in source and binary forms, with or without
%modification, are permitted provided that the following conditions are met :
%*Redistributions of source code must retain the above copyright
%notice, this list of conditions and the following disclaimer.
%*Redistributions in binary form must reproduce the above copyright
%notice, this list of conditions and the following disclaimer in the
%documentation and / or other materials provided with the distribution.
%* Neither the name of the <organization> nor the
%names of its contributors may be used to endorse or promote products
%derived from this software without specific prior written permission.
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
%ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
%WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
%DISCLAIMED.IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
%DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
%(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
%LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
%ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
%(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
%SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%

%k-step recover puncture pattern search
function [P, Rates, Pairs] = GroupAndSort2(H, MomRate, shag)
[m, n] = size(H);
G = cell(1,1);
R = cell(1,1);
Gzero = [];
Rzero = [];
RowAndCol = [];
%_____Initialization___________________________________________________
k = 1;
Ginf = 1:n; % all variable nodes
Lambda = cell(1,n); %nonzero positions in cols
for i = 1:n
    Lambda{i} = find(H(:,i));
end    
Rinf = 1:m; %all check nodes
Gamma = cell(1,m); %nonzero positions in rows
for i = 1:m
    Gamma{i} = find(H(i,:));
end
S = zeros(1, n);
count = 1;
Pairs = zeros(2,1);
%__Cycle__________________________________________________________________
    while ~isempty(Ginf)
        G(k)={[]};
        R(k)={[]};       
        while ~isempty(Rinf)
            %__STEP_1_STEP_2________
            Qinf = cell(1,length(Rinf));
            Qpower = zeros(1,length(Rinf));
            Qnumber = zeros(1,length(Rinf));
            for i = 1:length(Rinf)
                thisRow = Rinf(i);
                Qinf{i} = intersect(Gamma{thisRow},Ginf);
                Qpower(i) = length(Qinf{i});
                Qnumber(i) = thisRow;
            end
            minRWeff = min(Qpower);
            OmegaSlots = find(Qpower==minRWeff);
            Omega = Qnumber(OmegaSlots);
            
            %__STEP_3_STEP_3.1_STEP3.2____
            rows4Cinf = Qinf(OmegaSlots);
            cols4Cinf = [];
            for i = 1:length(rows4Cinf)
                cols4Cinf = [cols4Cinf rows4Cinf{i}];
            end
            cols4Cinf = setdiff(cols4Cinf,n+1);
            
            Cinf = cell(1,length(cols4Cinf));
            Cpower = zeros(1,length(cols4Cinf));
            for i = 1:length(cols4Cinf)
                Cinf{i} = intersect(Lambda{cols4Cinf(i)},Rinf);
                Cpower(i) = length(Cinf{i});
            end
            minCWeff = min(Cpower);
            minCols = cols4Cinf(find(Cpower==minCWeff));
            Omega2 = [];
            SmallC2 = [];
            for i = 1:length(rows4Cinf)
                cols2pair = intersect(rows4Cinf{i},minCols);
                if ~isempty(cols2pair)
                    Omega2 = [Omega2 Omega(i)];
                    point = length(cols2pair);
                    SmallC2 = [SmallC2 cols2pair(point)];
                end                
            end
            
            %__STER_3.3____________
            W = zeros(1,length(Omega2));
            for i = 1:length(Omega2)
               varNodes = Gamma{Omega2(i)};
               for j = 1:length(varNodes)
                   W(i) = W(i) + S(varNodes(j));
               end
            end
            minOfW = find(W==min(W));
            point = randi(length(minOfW),1,1);
            %point = randi(length(minOfW),1,1);
            ChosenRow = Omega2(point);
            ChosenCol = SmallC2(point);
            RowAndCol = [RowAndCol [ChosenRow;ChosenCol]];
            
            %__STEP_4.0__UPDATING___
            G{k}    = [G{k} ChosenCol];
            R{k}    = [R{k} ChosenRow];
            Pairs(1,count) = ChosenCol;
            Pairs(2, count) = ChosenRow;
            count = count + 1;
            
            exluded = setdiff(Qinf{find(Rinf==ChosenRow)},ChosenCol);
            
            Gzero   = [Gzero exluded];
            thisOfC = find(cols4Cinf==ChosenCol);
            Rzero   = [Rzero setdiff(Cinf{thisOfC},ChosenRow)'];
            
            Ginf    = setdiff(Ginf,Qinf{find(Rinf==ChosenRow)});
            Rinf    = setdiff(Rinf,Cinf{thisOfC});
            
            temp = exluded;
            for i = 1:length(temp)
                S(temp(i)) = 1;
            end
            
            temp = setdiff(Gamma{ChosenRow},ChosenCol);
            S(ChosenCol)=0;
            for i = 1:length(temp)
                S(ChosenCol) = S(ChosenCol) + S(temp(i));
            end
            %__STEP_5.0__CHECK____
            if isempty(Ginf)
                break
            end
        end
        Rzero = setdiff(Rzero,m+1);
        Gzero = setdiff(Gzero,n+1);
        
        %__STEP_5.1_STEP5.2_ADD_CHECK____        
        if (~isempty(Ginf))&&(isempty(Rinf))
            testWReff = zeros(1,length(Rzero));
            for i = 1:length(Rzero)
                testWReff(i) = RWeff(Rzero(i),Ginf,H);
            end                
            Rinf = Rzero(find(testWReff));                
        end
        
        %__STEP_6.0___
        k = k + 1;
               
    end

    
    %__________________Sorting Part of Algorithm_______________________________
    %__________________________________________________________________________

    %_________Declarations and Definitions_____________________________________
    MaxIter = k-1;
    sum=0;
    for i = 1:MaxIter
        sum = sum + length(G{i})/n;
    end
%     MomRate = (n-m)/n;
    MaxRate = MomRate/(1-sum);
    Rates = MomRate:shag:MaxRate;
    if Rates(end)<MaxRate
        Rates = [Rates MaxRate];
    end
    M = length(Rates)-1;
   
    %_________Initialization___________________________________________________
    R = 1:m;
    j = 2;
    k = 1;
    P = {[]};
    %_________Cycle____________________________________________________________
    while j<=M+1
        P{j}=P{j-1};
        npj = fix(n*(Rates(j) - MomRate)/Rates(j));
        dnpj = npj - length(P{j});
        while dnpj~=0
            CW = zeros(1,length(G{k}));
            for i = 1:length(G{k})
                CW(i) = CWeff(G{k}(i), R, H);
            end
            porog = max(CW);
            C = G{k}(find(CW==porog));
            if length(C)>1
                degC = zeros(1,length(C));
                for i = 1:length(C)
                    degC(i) = CWeff(C(i),1:m,H);
                end
                porog = min(degC);
                chC = find(degC==porog);
                c2 = C(randi(length(chC)));
                %c2 =L*Z+1;
                P{j} = [P{j} c2];
                G{k} = setdiff(G{k},c2);
                R = setdiff(R, Lambda{c2});
            elseif length(C)==1
                P{j} = [P{j} C];
                G{k} = setdiff(G{k}, C);
                R = setdiff(R, Lambda{C});
            end

            dnpj = dnpj - 1;
            if isempty(G{k})
                k = k+1;
                R = 1:m;
            end
        end
        j = j + 1;
    end

end


function M_ = SeamlessParametrization(M,l0,thetaT,epsilon)    
    l = 2*log(l0);  
    
    C = initConstraints(M);
    ThetaT = [thetaT.vertexAngles;thetaT.loopAngles];
    while true
        [M_, l_, D, C] = DiffMakeDelaunay(M, l, C);
        
        [alpha, grad_alpha] = ComputeAnglesAndGradient(M_, l_);
        assignin('base','alpha',alpha);
        assignin('base','grad_alpha',grad_alpha);
        F = C*alpha-ThetaT;
        gradF = C * grad_alpha * D;
        assignin('base','F',F);
        assignin("base","gradF",gradF);
        L = gradF * gradF';
        assignin('base',"L",L);
        lambda = 1e-6;
        L_reg = L + lambda * eye(size(L));
        mu = L_reg \ (-F);
        % mu = L\-F;
        assignin('base','mu',mu);
        d = gradF' * mu;
        assignin('base','d',d);
        beta = LineSearch(M_, C, ThetaT, l, d, F);
        assignin('base','beta',beta);
        assignin('base','l_old',exp(l/2));
        l = l + beta * d;
        assignin('base','l_new',exp(l/2));
        if norm(C*alpha - ThetaT,'inf') < epsilon
            disp('yes');
            break;
        end
        disp(norm(C*alpha - ThetaT,'inf'));
    end
end

function [M_,l_,D,C] = DiffMakeDelaunay(M,l,C)
    M_ = M; l_ = l;
    E = M.edges; NE = size(E,1);
    D = eye(NE);
    
    deque = 1:NE;
    while true
        ei = deque(1);
        e = E(ei,:);
        assignin('base','ei',ei);
        assignin('base','e',e);
        if NonDelaunay(M_,ei,e)
            [M_,l_,e, eiabcd] = PtolemyFlip(M_,l_,ei,e);
            
            D = DiffPtolemy(M_,l_,ei,e) * D;
            
            C = UpdateConstraints(M_,l_,C,e);
            
            deque = [deque, ei, eiabcd];
        end
        deque(1) = [];

        if isempty(deque)
            break;
        end
    end
end

function C = initConstraints(M)
    NV = size(M.vertices, 1); 
    NF = size(M.faces, 1);    
    g = M.genus;     

    C = zeros(NV + 2*g - 1, 3*NF);

    for vi = 1:NV - 1
        [fi, i] = find(M.faces == vi);
        C(vi, 3*(fi-1) + i) = 1; 
    end

    for j = 1:length(M.loops)
        loop = M.loops{j}; 
        [d_j, k] = computeLoopSigns(loop, M.Fe, M.edges, M.faces); 
        for m = 1:length(loop) 
            fi = loop(m); 
            C(NV + j - 1, 3*(fi - 1) + k(m)) = d_j(m); 
        end
    end
end

function [d_j, k] = computeLoopSigns(loop, Fe, E, F)
    d_j = zeros(length(loop), 1);
    k = zeros(length(loop), 1);
    for m = 1:length(loop)
        fi = loop(m); 
        fip = loop(mod(m - 2, length(loop)) + 1);
        fin = loop(mod(m, length(loop)) + 1);
        assignin('base','cLS1',[m length(loop) fip fi fin]);
        assignin('base','cLS2',Fe([fip fi fin],:));
        ei1 = intersect(Fe(fip, :), Fe(fi, :));
        ei2 = intersect(Fe(fi, :), Fe(fin, :));
        d_j(m) = determineSign(fi, ei1, ei2, Fe);

        vi = intersect(E(ei1,:), E(ei2,:));
        k(m) = find(F(fi,:) == vi);
    end
end

function sign = determineSign(fi, ei1, ei2, Fe)
    assignin('base','fi',fi);
    assignin('base','ei1_dS',ei1);
    assignin('base','ei2_dS',ei2);
    assignin('base','Fe',Fe);
    i = find(Fe(fi,:) == ei1);
    if Fe(fi,mod(i,3)+1) == ei2
        sign = 1;
    else
        sign = -1;
    end
end

function isNonDelaunay = NonDelaunay(M,ei,e)
    Fe = M.Fe;
    l = M.edgeLengths;
    fidxes = M.TR.edgeAttachments(e);
    fi1 = fidxes{:}(1); fi2 = fidxes{:}(2);
    
    abe = Fe(fi1,:); cde = Fe(fi2,:);
    ab = abe(~ismember(abe,ei)); cd = cde(~ismember(cde,ei));
    a = ab(1); b = ab(2); c = cd(1); d = cd(2);
    cos1 = (exp(l(a))+exp(l(b))-exp(l(ei)))/(2*exp((l(a)+l(b))/2));
    cos2 = (exp(l(c))+exp(l(d))-exp(l(ei)))/(2*exp((l(c)+l(d))/2));
    
    if cos1 + cos2 >= 0
        isNonDelaunay = false;
    else
        isNonDelaunay = true;
    end
end

function [M_,l_,e_,eiabcd] = PtolemyFlip(M_,l_,ei,e)
    F = M_.faces; F_ = F;
    E = M_.edges; E_ = E;
    Fe = M_.Fe; Fe_ = Fe;
    fidxes = edgeAttachments(M_.TR, e);
    fi1 = fidxes{:}(1); fi2 = fidxes{:}(2);
    assignin('base','fidxes_PF',fidxes);
    
    v1 = F(fi1, ~ismember(F(fi1,:), e));
    v2 = F(fi2, ~ismember(F(fi2,:), e));
    i11 = find(F(fi1,:) == v1); F_(fi1,mod(i11+1,3)+1) = v2; v31 = F(fi1,mod(i11,3)+1);
    i22 = find(F(fi2,:) == v2); F_(fi2,mod(i22+1,3)+1) = v1; v32 = F(fi2,mod(i22,3)+1);
    
    e_ = sort([v1 v2]);
    E_(ei,:) = e_;
    
    a = sort([v1 v31]); eia = find(all(E_ == a,2));
    b = sort([v31 v2]); eib = find(all(E_ == b,2));
    c = sort([v2 v32]); eic = find(all(E_ == c,2));
    d = sort([v32 v1]); eid = find(all(E_ == d,2));
    ia1 = find(Fe(fi1,:) == eia); Fe_(fi1,mod(ia1,3)+1) = eib; Fe_(fi1,mod(ia1+1,3)+1) = ei;
    ic2 = find(Fe(fi2,:) == eic); Fe_(fi2,mod(ic2,3)+1) = eid; Fe_(fi2,mod(ic2+1,3)+1) = ei;
    eiabcd = [eia eib eic eid];
    save('data.mat','v1','v2','v31','v32','fi1','fi2','a','b','c','d','eia','ia1');

    l_(ei) = 2*log((exp((l_(eia)+l_(eic))/2)+exp((l_(eib)+l_(eid))/2)))-l_(ei);    
    
    assignin('base','M_pre',M_);
    for j = 1:length(M_.loops)
        loop = M_.loops{j};
        ll = length(loop);
        if any(loop == fi1) && ~any(loop == fi2)
            m1 = find(loop == fi1);
            if any(Fe(loop(mod(m1, ll)+1),:) == eia)
                loop = [loop(1:m1-1) fi2 fi1 loop(m1+1:end)];
            elseif any(Fe(loop(mod(m1, ll)+1),:) == eid)
                loop = [loop(1:m1-1) fi1 fi2 loop(m1+1:end)];
            end
        elseif ~any(loop == fi1) && any(loop == fi2)
            m2 = find(loop == fi2);
            if any(Fe(loop(mod(m2, ll)+1),:) == eib)
                loop = [loop(1:m2-1) fi2 fi1 loop(m2+1:end)];
            elseif any(Fe(loop(mod(m2, ll)+1),:) == eic)
                loop = [loop(1:m2-1) fi1 fi2 loop(m2+1:end)];
            end
        elseif any(loop == fi1) && any(loop == fi2)
            m1 = find(loop == fi1); m2 = find(loop == fi2);
            mm = min(m1,m2); mM = max(m1,m2);
            eidxes = [Fe(loop(mod(mm-2,ll)+1),:) Fe(loop(mod(mM,ll)+1),:)];
            if any(eidxes == eia) && any(eidxes == eib)
                loop(loop == fi2) = [];
            elseif any(eidxes == eib) && any(eidxes == eid)
                loop([m1 m2]) = loop([m2 m1]);
            elseif any(eidxes == eic) && any(eidxes == eid)
                loop(loop == fi1) = [];
            end
        end
        M_.loops{j} = loop;
    end

    M_.faces = F_;
    M_.TR = triangulation(F_,M_.vertices);
    M_.edges = E_;
    M_.Fe = Fe_;
    M_.edgeLengths = l_;

    assignin('base','M_',M_);
end

function dD = DiffPtolemy(M_,l_,ei,e)
    Fe = M_.Fe;
    fidxes = M_.TR.edgeAttachments(e);
    assignin('base',"fidxes_DP",fidxes);
    fi1 = fidxes{:}(1); fi2 = fidxes{:}(2);
    assignin('base','fi1_DP',fi1);
    assignin('base','fi2_DP',fi2);
    
    i1 = find(Fe(fi1,:) == ei);
    a = Fe(fi1,mod(i1,3)+1); b = Fe(fi1,mod(i1-2,3)+1);
    i2 = find(Fe(fi2,:) == ei);
    c = Fe(fi2,mod(i2,3)+1); d = Fe(fi2,mod(i2-2,3)+1);
    
    t = exp((l_(a)+l_(c)-l_(b)-l_(d))/2);
    
    dD = eye(length(l_));
    assignin('base','De',[ei a b c d]);
    dD(ei,[ei a b c d]) = [-2 2*t/(1+t) 2/(1+t) 2*t/(1+t) 2/(1+t)];
end

function C = UpdateConstraints(M_,l_,C,e)
    C = initConstraints(M_);
end

function [alpha, grad_alpha] = ComputeAnglesAndGradient(M,l)
    Fe = M.Fe;
    NF = size(Fe, 1);
    NE = size(l, 1);
    alpha = zeros(3*NF, 1);
    grad_alpha = zeros(3*NF, NE);  
    
    for fi = 1:NF
        ei1 = Fe(fi, 1); ei2 = Fe(fi, 2); ei3 = Fe(fi, 3);
        l1 = exp(l(ei1)/2); l2 = exp(l(ei2)/2); l3 = exp(l(ei3)/2);
        
        u1 = (l2^2 + l3^2 - l1^2) / (2*l2*l3);
        u2 = (l1^2 + l3^2 - l2^2) / (2*l1*l3);
        u3 = (l1^2 + l2^2 - l3^2) / (2*l1*l2);

        alpha(3*(fi-1)+1) = acos(u1);
        alpha(3*(fi-1)+2) = acos(u2);
        alpha(3*(fi-1)+3) = acos(u3);
        
        grad_alpha(3*(fi-1)+1,Fe(fi,1)) = l1/(l2*l3*sqrt(1-u1^2));
        grad_alpha(3*(fi-1)+1,Fe(fi,2)) = -(l1^2+l2^2-l3^2)/(2*l2^2*l3*sqrt(1-u1^2));
        grad_alpha(3*(fi-1)+1,Fe(fi,3)) = -(l1^2+l3^2-l2^2)/(2*l2*l3^2*sqrt(1-u1^2));
  
        grad_alpha(3*(fi-1)+2,Fe(fi,1)) = -(l1^2+l2^2-l3^2)/(2*l1^2*l3*sqrt(1-u2^2));
        grad_alpha(3*(fi-1)+2,Fe(fi,2)) = l2/(l1*l3*sqrt(1-u2^2));
        grad_alpha(3*(fi-1)+2,Fe(fi,3)) = -(l3^2+l2^2-l1^2)/(2*l1*l3^2*sqrt(1-u2^2));
            
        grad_alpha(3*(fi-1)+3,Fe(fi,1)) = -(l1^2-l2^2+l3^2)/(2*l1^2*l2*sqrt(1-u3^2));
        grad_alpha(3*(fi-1)+3,Fe(fi,2)) = -(l2^2-l1^2+l3^2)/(2*l1*l2^2*sqrt(1-u3^2));
        grad_alpha(3*(fi-1)+3,Fe(fi,3)) = l3/(l1*l2*sqrt(1-u3^2));
    end
end

function beta = LineSearch(M, C, ThetaT, l, d, F)
    beta = 1;
    max_iter = 500; 

    for iter = 1:max_iter
        l_new = l + beta * d;
        [alpha_new, ~] = ComputeAnglesAndGradient(M, l_new);
        F_new = C * alpha_new - ThetaT;
        
        if norm(F_new) <= norm(F) && dot(F,F_new) >= 0
            return;
        end
        
        beta = beta * 0.9;
    end
    beta = 1.0;
end

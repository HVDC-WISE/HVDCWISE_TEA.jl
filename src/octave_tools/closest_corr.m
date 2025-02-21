function COR1 = closest_corr(COR)
    [V D] = eig(COR,'nobalance');
    Dp = D;

    for jd = 1:size(D,1)
        if D(jd,jd) <= 1e-2
            Dp(jd,jd) = 1e-2;
        end
    end
    T = zeros(size(D));

    eigenv_p = diag(Dp);

    for i = 1:size(T,1)
        T(i,i)= 1/((V(i,:).^2)*eigenv_p); 
    %      T(i,i) = dummy;
    end

    Bp = V * sqrt(Dp);

    Bf = sqrt(T)*Bp;

    COR1 = Bf*Bf';
    if COR(1,1)< 0.1
        COR1 = COR1+eye(size(COR1));    
    end
end
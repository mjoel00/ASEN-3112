function [] = ploteigenvector(L,ev,ne,nsub,scale)
nv = ne*nsub+1;
    Le = L/ne;
    dx = Le/nsub;
    x = zeros(nv,1);
    v = x;
    k=1;
    for i = 1:ne
        xi = Le*(i-1);
        vi = ev(2*i-1);
        qi = ev(2*i);
        vj = ev(2*i+1);
        qj = ev(2*i+2);
        for n = 1:nsub
            xk = xi+dx*n;
            zeta = (2*n-nsub)/nsub;
            vk = scale*(.125*(4*(vi+vj)+2*(vi-vj)*((zeta^2)-3)*zeta+...
                Le*((zeta^2)-1)*(qj-qi+(qi+qj)*zeta)));
            k = k+1;
            x(k) = xk;
            v(k) = vk;
        end
    end

        figure
        plot(x,v,'Linewidth',2)
        title('Mode Shape')
        grid on
  
   
    %saveas(gcf,'4ElementModalBehavior.jpg','jpg')
end
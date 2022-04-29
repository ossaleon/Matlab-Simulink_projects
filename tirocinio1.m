clear
clc
nu = 5;
lambda_1 = 1;
lambda_2 = ;
B = [1; 1];
A = [lambda_1 0; 0 lambda_2];

%grafici(A, B, nu, [sqrt(2)/2; sqrt(2)/2])
grafici(A, B, nu, [1; -1])
figure(4)
min_energia(A, B, nu)
grid on

%calcola u ottimo dati A, B, nu, x segnato
function u = controllo_ottimo(A, B, nu, x_segn)
    P = B;
    %se nu > 1 calcola P in maniera iterativa, altrimenti P=B
    %da rivedere preallocazione matrice
    if nu > 1
        temp = B;
        for i = 1:nu-1
            temp = A*temp;
         P = [P, temp];
        end
    end
    
    G = P*P'; %calcolo G
    if det(G) == 0 %se G non è invertibile risolve il sistema G*beta = x
        %scegliendo un beta con variabile libera beta1 = 1
        syms beta2
        beta1 = 1;
        beta = [beta1; beta2];
        eq = G*beta;
        eqs = [eq(1,:) == x_segn(1,:); eq(2,:) == x_segn(2,:)];
        beta = [beta1, solve(eqs)];
        if size(beta) == 1 %se il sistema non è risolvibile poichè x_segn non
            %è raggiungibile, il controllo ottimo non esiste e ritorna un vettore di Nan
            [~, colnum] = size(B);
            u = NaN(nu*colnum, 1);
        end    
    else
        beta = G\x_segn; %risoluzione sistema lineare inv(G)*x_segn = beta
        u = flip(P'*beta); %uso flip per calcolare il vettore w e poi averlo al contrario
    end
end

%realizza i grafici sovrapposti (al variare di nu) alpha-J
function min_energia(A, B, nuMax)
    %crea un vettore di 100 angoli in radianti compresi tra 0 e 2*pi equidistanziati
    alpha = linspace(0, 2*pi); 
    
    %equazioni parametriche della circonferenza unitaria
    x1 = cos(alpha);
    x2 = sin(alpha);
    x = [x1', x2']'; %matrice che ha per colonne x = [x1; x2], ovvero i punti della circonferenza unitaria
    
    [~, numCols] = size(x);
    %realizza un grafico alpha-J per ogni nu in [0, nuMax]
    for nu = 1:nuMax
        J = [];
        for i = 1:numCols %per ogni punto della circonferenza preso in considerazione
            u = controllo_ottimo(A, B, nu, x(:,i)); %calcolo u ottimo
            temp = u'*u; %calcolo Jnu(u)
            J = [J, temp]; %vettore che contiene Jnu(u) al variare di u fissato un nu
        end
        %grafico ascissa alpha - ordinata J
        plot(alpha, J,'DisplayName', sprintf('nu = %d', nu))
        hold on %sovrappone i grafici Jnu(u) al variare di nu
    end
    hold off
    legend
    title('Andamento di J al variare di u* parametrizzato con nu')
    xlabel('alpha')
    ylabel('J')
end

%realizza i tre grafici nu-u*, nu-x*1;x*2, x*1-x*2
function grafici(A, B, nuMax, x_star)
    %calcolo u ottimo dati A, B, nuMax  
    temp = controllo_ottimo(A, B, nuMax, x_star);
    %creo un vettore riga composto da u e 0 (u(nu)=0)
    u_star = [temp', 0];
    
    %grafico nu ascissa - u* ordinata
    figure(1)
    plot([0:nuMax], u_star, '--ob')
    
    title('nu - u*')
    xlabel('nu')
    ylabel('u*')
    
    %storico di x nell'intervallo [0, nu]
    x_star = [0; 0]; %x(0)=0
    for nu = 1:nuMax
        temp = A*x_star(:,nu)+B*u_star(nu);
        x_star = [x_star, temp]; %ogni colonna corrisponde a x(nu-1). vedere preallocazione
    end
    
    %grafico nu ascissa - x*1, x*2 ordinata
    figure(2)
    plot([0:nuMax], x_star(1,:), '--om')
    hold on
    plot([0:nuMax], x_star(2,:), '--og')
    legend('x*_1', 'x*_2')
    title('nu - x*_1, x*_2')
    xlabel('nu')
    ylabel('x*_1, x*_2')
    hold off
    
    %grafico x*1 ascissa - x*2 ordinata
    figure(3)
    plot(x_star(1,:), x_star(2,:), '--or')
    title('x*_1 - x*_2')
    xlabel('x*_1')
    ylabel('x*_2')
    
end


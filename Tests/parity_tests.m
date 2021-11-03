%Tests for calls
call = 1; 
prices_ui = barrier_prices_cm(S_c, K, B, "uo", call); 
prices_uo = barrier_prices_cm(S_c, K, B, "ui", call); 
prices_di = barrier_prices_cm(S_c, K, B, "di", call); 
prices_do = barrier_prices_cm(S_c, K, B, "do", call); 

up = prices_ui + prices_uo;
down = prices_di + prices_do;

vannilla = vanilla_prices(S_c, K, 1);

allColEqual(up);
isalmost(up, down);
isalmost(up(:, 1), vannilla);



%Tests for puts
call = -1;
prices_ui = barrier_prices_cm(S_c, K, B, "uo", call); 
prices_uo = barrier_prices_cm(S_c, K, B, "ui", call); 
prices_di = barrier_prices_cm(S_c, K, B, "di", call); 
prices_do = barrier_prices_cm(S_c, K, B, "do", call); 

up = prices_ui + prices_uo;
down = prices_di + prices_do;

vannilla = vanilla_prices(S_c, K, call);

allColEqual(up);
isalmost(up, down);
isalmost(up(:, 1), vannilla);



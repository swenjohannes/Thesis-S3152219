function tenor_vector_out = str2tenor(Ymwd,base)

if nargin < 2
    base = 0;
end

Ymwd = cellstr(Ymwd);
%tenor_vector_out = [];
tenor_vector_out = zeros(length(Ymwd),6);
if length(base) == 1 % scalar
    base = base*ones(size(Ymwd));
elseif length(base) ~= length(Ymwd)
    error('Input vectors must agree in size.');
end
for j = 1:length(Ymwd)
    ymwd = char(Ymwd(j));
%    [b,e] = regexpi(ymwd,{'-{0,1}\d*y', '-{0,1}\d*m', '-{0,1}\d*w', '-{0,1}\d*d', '-{0,1}\d*[obt]'});
    [b,e] = regexpi(ymwd,{'-{0,1}\d*y', '-{0,1}\d*m', '-{0,1}\d*w', '-{0,1}\d*d{1}', '-{0,1}\d*[obt]', '-{0,1}\d*d{2}','-{0,1}\d*s{1}.{0,1}n{1}','-{0,1}\d*c{1}','-{0,1}sw'});
    % -{0} = ze
    tenor_vector = zeros(1,5);
    
    for i = 1:5
        if ~isempty(ymwd(b{i}:e{i}-1))
            tenor_vector(i) = str2double(ymwd(b{i}:e{i}-1));
        end
    end
    if ~isempty(b{5}) && isempty(ymwd(b{5}:e{5}-1)) % Treat 'O/N' and 'OD'case
        tenor_vector(5) = tenor_vector(5) + 1;
    end
    if ~isempty(b{6}) % Treat 'DD'case
        tenor_vector(5) = tenor_vector(5) + 1;
    end
    if ~isempty(b{7}) % Treat 'S/N'case
        tenor_vector(5) = tenor_vector(5) + 1;
    end
    if ~isempty(b{8}) % Treat 'C'case
        tenor_vector(5) = str2double(ymwd(b{8}:e{8}-1));
    end
    if ~isempty(b{9}) % Treat 'SW' case
        tenor_vector(3) = 1;
    end
    tenor_vector = [tenor_vector(1:2), tenor_vector(3:5)*[7;1;1], [0 0 0]] + datevec(base(j));
    %tenor_vector_out = [tenor_vector_out; tenor_vector];
    tenor_vector_out(j,:) = tenor_vector;
end


% http://www.finansraadet.dk/english/menu/CIBOR/Tomorrow+Next+T+N/
% T/N Tomorrow Next interest rates
% 
% Rules for calculation of the Tomorrow/Next interest rate
% 
% The Tomorrow/Next interest rate
% The tomorrow/next interest rate (T/N) is a reference rate for unsecured money market lending (deposit lending) in Danish kroner (DKK). The transaction starts on the first banking day following the trading day and matures on the second banking day following the trading day.
% 
% Technical specifications 
% The calculation of the T/N interest rate is based on daily reports from a number of banks:
% 
%     * Each declaring bank reports the amount of unsecured T/N inter-bank lending/deposits denominated in Danish Kroner and the weighted average interest rate for these lending transactions.
%     * Data is reported on all Danish banking days. Turnover outside Danish banking days shall not be reported.
%     * Lending transactions comprise transactions involving Danish and for-eign banks traded directly or through a broker. Transactions between branches and subsidiaries of a credit institution and transactions with Danmarks Nationalbank are not included in the reported transactions.
%     * Data is reported with a time lag of one day, e.g. Monday's lending transactions are reported on Tuesday at the latest. Data is reported to Danmarks Nationalbank at 10.00 a.m. CET at the latest covering in-formation on lending transactions from the previous banking day. The daily information includes lending transactions traded between 8.00 a.m. and 5.00 p.m. CET. The sum of lending transactions shall be re-ported by banks in millions of Danish kroner and the individual interest rate is computed as a weighted average with 2 decimals.
%     * The T/N interest rate is computed according to the interest rate convention on the Danish money market, i.e. ACT/360.
%     * Based on the reported individual data, Danmarks Nationalbank calculates an average interest rate in which the interest rate of each bank is weighted according to the bank's share of total reported lending trans-actions.
% 
% Reporting institutions
% The current list of T/N reporting institutions:
% 
% Publication
% The T/N interest rate is published on the website of Danmarks Nationalbank at 12.00 on the day of reporting. The interest rate is published with 4 deci-mals.
% 
% The T/N interest rate is also shown on the Danish Bankers Associations (Finansrådet) homepage.
% 
% The T/N Committee
% The Danish Bankers Association (Finansrådet) administrates the T/N, whereas Danmarks Nationalbank performs the daily T/N computation.
% 
% In connection with the T/N reference rate, a Committee has been established under the Danish Bankers Association in order to address issues concerning the T/N interest rate, including determining the list of reporting banks. Each reporting bank is a member of the Committee. Danmarks Nationalbank participates as an observer.

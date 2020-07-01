t_rect=zeros(50,1);
t_sparse=zeros(50,1);

for L=1:100
% t = cputime;
% PCF_normalisation_rectangular(L,L);
% e = cputime-t;
% t_rect(L)=e;

t = cputime;
PCF_normalisation_rectangular_sparse_JO(L,L);
e = cputime-t;
t_sparse(L)=e;
end

plot(t_sparse);
% hold on
% plot(t_rect);

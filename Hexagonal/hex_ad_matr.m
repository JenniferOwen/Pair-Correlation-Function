%Thsi code generate an irregulare lattice starting from a regular square
%lattice, perturbating the coordinates at random and then considering the
%Varanoi Tasselation. The adjacency matrix is computed.
%Enrico Gavagnin 28/07/2017

function [NORM]=hex_ad_matr(Lx,Ly)

close all
    %Label all sites
    all_numb=0;
    
    %Label occupied sites
    occ_numb=0;
    
    %visualisation
    PLOT_ON=1;
    
    MAP=[1 1 1;.3 .3 .3];
    colormap(MAP);
    caxis([0 1]);
    %Specify the density of occupued cells.
    density=0;
    %Dimension of the original regular lattice (square)
    
    
    x=[];
    y=[];
    for i=2:Lx-1
         for j=2:Ly-1
%                      x=[x;i+1.3*(rand-0.5)];
%                      y=[y;j+1.3*(rand-0.5)];
            if mod(i,2)==0
               y=[y;j+0.5];
            else
               y=[y;j];
            end
            x=[x;i];
        end
    end
    %coord = gallery('uniformdata',[10 2],5);
    coord=[x,y];
    
    %Varanoi Tasselation
    [v,c] = voronoin(coord);
    
    %Initialisation of some indeces (total number of triangles, occupied
    %triangles and counter of total triangles
    index=[];
    index_occ=[];
    tot_trian=0;
    
    %% Visualisation
    
    for i = 1:length(c)
        if all(c{i}~=1)% If at least one of the indices is 1,
            % then it is an open region and we can't
            % patch that.
            
            %We esclude also some triangles with highly irregular shape, in
            %particular if their vertices are outide the original regular lattice
            if all(all(v(c{i},:)>1))&& (all(all((v(c{i},:)<Lx))| all((v(c{i},:)<Ly))))  
                tot_trian=tot_trian+1;
                index=[index;i];
                
                %Occupy at random some of the sites
                if rand<density;
                    %Colour of the occupied triangles
                    col=2;
                    index_occ=[index_occ;tot_trian];
                    if PLOT_ON
                        %display the i-th triangle
                        patch(v(c{i},1),v(c{i},2),col);
                        if occ_numb
                            %Write the number of the triangle in the centre of it
                            text(mean(v(c{i},1)),mean(v(c{i},2)),sprintf('%d',size(index_occ,1)))
                        end
                    end
                else
                    %Colour of the empty triangles
                    col=0;
                    if PLOT_ON
                        %display the i-th triangle
                        patch(v(c{i},1),v(c{i},2),col,'LineWidth',2);
                    end
                end
                
                if all_numb && PLOT_ON
                    %Write the number of the triangle in the centre of it
                    text(mean(v(c{i},1)),mean(v(c{i},2)),sprintf('%d',tot_trian),'FontSize',24)
                end
            end
        end
    end
    
    %% Adjacency Matrix
    
    % Initialisation
    A=zeros(size(index,1));
    
    
    
    for tri1=1:size(index,1)-1
        for tri2=tri1+1:size(index,1)
            
            % Determine if the two triangles tri1 and tri2 are neighbours: if they
            % have at least two vertices in common
            if size(intersect(c{index(tri1)},c{index(tri2)}),2)>1
                
                %In case they are neighbours insert a connection
                A(tri1,tri2)=1;
            end
        end
    end
    storeA=NaN(size(A,1),size(A,2),2*min(Lx,Ly));
    
    %Store the identity matrix
    storeA(:,:,1)=eye(size(A));
    
    
    %A is symmetric since if a is connected to b, b is connected to a.
    storeA(:,:,2)=A+A';
    i=2;
    
    %% Computation of the distance matrix and the powers of the adjacency matrix
    
    
    %Initialisation
    Dist=eye(size(A));
    storeA=NaN(size(A,1),size(A,2),2*min(Lx,Ly));
    
    %Store the identity matrix
    storeA(:,:,1)=eye(size(A));
    
    %A is symmetric since if a is connected to b, b is connected to a.
    storeA(:,:,2)=A+A';
    i=2;
    
    %Loop that compute the powers of the A. matrix inductively and
    %automatically update the matrix of the distances
    
    while min(min(storeA(:,:,i-1)))==0 && i<=2*min(Lx,Ly)
        
        %Update the matrix that store the pairwise distances
        Dist=Dist+(i-1)*double(storeA(:,:,i)>0).*double(Dist==0);
        i=i+1;
        
        %Compute the i+1-th power of A
        storeA(:,:,i)=storeA(:,:,i-1)*storeA(:,:,2);
    end
    %Correct the matrix of distances
    Dist=Dist-eye(size(A));
    
    %store the maximum path length
    max_path=i;
   NORM=[]; 
for i=1:max(max(Dist))
NORM=[NORM, nnz(Dist==i)];
end

end


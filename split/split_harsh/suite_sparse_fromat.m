
n=4;

A=(triu(rand(n,n)));
A=A+A';


%% extract the indices

At=(triu(A)); %%% find the vector with non zero elements

valuetot=nnz(At);

value_vec=zeros(valuetot,1);
column_index=zeros(valuetot,1);
row_index=zeros(n,1);
sum_element=0;

for k=1:n
    
    [i,j,s]=find(At(:,k));
    
    
    column_index(sum_element+1:sum_element+length(s))=i;
    value_vec(sum_element+1:sum_element+length(s))=s;
    
    sum_element=sum_element + length(s);
    row_index(k)=sum_element;   
    
    
end

%% Now save in the suite sparse format i.e. row index starting from zero

row_index=row_index-1;
column_index=[0;column_index];


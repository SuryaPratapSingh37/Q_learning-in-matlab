%% Q-learning Algorithm to solve a Grid-world
% *Authors:*Surya Pratap Singh, Yug Ajmera
% 

%% Size of square shaped grid-world
n = 10;
%% Assigning the rewards of -1 for normal states and -10 for obstacles

maze = -1*ones(n,n);

maze(1,4)=-10;
maze(3,7)=-10;
maze(3,8)=-10;
maze(3,9)=-10;
maze(4,3)=-10;
maze(4,7)=-10;
maze(5,3)=-10;
maze(5,5)=-10;
maze(5,8)=-10;
maze(6,3)=-10;
maze(6,4)=-10;
maze(6,5)=-10;
maze(7,7)=-10;
maze(7,8)=-10;
maze(8,4)=-10;
maze(9,4)=-10;
maze(10,4)=-10;
maze(10,7)=-10;


%% 
% * Starting position reward

maze(1,1) = -1;
%% 
% * Goal position reward

maze(n,n) = 100;

%% To display the grid-world

disp(maze)

n=length(maze);

figure
imagesc(maze)
colormap(white)

for i=1:n
    for j=1:n
        if maze(i,j)==-10
            rectangle('Position',[j-0.2 i-0.2 1.0 1.0],'FaceColor',[0 .5 .5]);%text(j,i,'X','HorizontalAlignment','center')
        end
    end
end
cx=0;
for i=1:n
    for j=1:n
        if (maze(i,j)==-1 && cx>0)
            rectangle('Position',[j-0.2 i-0.2 1.0 1.0],'FaceColor',[1 1 1]);%text(j,i,'X','HorizontalAlignment','center')
        end
        cx=cx+1;
    end
end
text(1+0.25,1+0.25,'Start','HorizontalAlignment','center')
text(n+0.25,n+0.25,'Goal','HorizontalAlignment','center')

axis off

Goal=n*n;
fprintf('Goal State is: %d',Goal)



%% Reward Matrix
% 
%%
reward=maze;

ab=zeros(n*n);


%% Q-Learning algorithm

%%
q = zeros(size(ab));%randn(size(ab));
gamma = 0.9;
alpha = 0.8;
maxItr = 200;
epsilon=0.9;
epsilon_discount = 0.9;
%% 

smatrix=zeros(n);
t=0;
for i=1:n
    for j=1:n
        smatrix(i,j)=t+1;
        t=smatrix(i,j);
    end
end

l=0;
select=0;
cum_matrix=[0 0];
epi_matrix=[0 0];
step_matrix=[0 0];
for i=1:maxItr
    epsilon = epsilon * epsilon_discount;
    cumulative = 0;
    % Starting from start position    
    cs=1;
    cs_row=1;
    cs_col=1;
    
    l=l+1;
    fprintf('Episode no.=%d\n',l);
    sd=0*ones(50);
    u=0;
    
    % Repeat until Goal state is reached
    while(1)
        u=u+1;
        %fprintf('Step no.=%d\n',u);
        
        sd(u,u)=cs;
    % possible actions for the chosen state
    if (cs_row==1)&&(cs_col==1)
        n_actions =[cs_col+1 n+1];
        
        %fprintf('Condition 1\n');
    end
    if (cs_row==1)&&(cs_col>1 && cs_col<n)
        n_actions =[cs_col-1 cs_col+1 cs_col+n];
        %fprintf('Condition 2\n');
    end
    if (cs_row==1)&&(cs_col==n)
        n_actions =[cs_col-1 2*n];
        %fprintf('Condition 3\n');
    end
    if (cs_row>1 && cs_row<n) && (cs_col==n)
        n_actions =[n*(cs_row-1) (n*cs_row-1) n*(cs_row+1)];
        %fprintf('Condition 4\n');
    end
    if (cs_row==n) && (cs_col==n)
        n_actions =[n*(n-1) (n*n-1)];
        %fprintf('Condition 5\n');
    end
    if (cs_row==n) && (cs_col<n && cs_col>1)
        n_actions =[n*(cs_row-1)+cs_col+1 n*(cs_row-1)+cs_col-1 n*(cs_row-1)+cs_col-n];
        %fprintf('Condition 6\n');
    end
    if (cs_row==n) && (cs_col==1)
        n_actions =[n*(cs_row-1)+2 n*(cs_row-1)+1-n];
        %fprintf('Condition 7\n');
    end
    if (cs_row<n && cs_row>1) && (cs_col==1)
        n_actions =[n*(cs_row-1)+2 n*(cs_row-1)+1-n n*(cs_row-1)+1+n];
        %fprintf('Condition 8\n');
    end
    if (cs_row<n && cs_row>1) && (cs_col<n && cs_col>1)
        n_actions =[n*(cs_row-1)+cs_col+1 n*(cs_row-1)+cs_col-1 n*(cs_row-1)+cs_col-n n*(cs_row-1)+cs_col+n];
        %fprintf('Condition 9\n');
    end
    
        
    % choose an action at random and set it as the next state
    z=rand();
    if (epsilon>z)
        select=1;
    end
    
    if (epsilon<=z)
        select=2;
    end
    
    %fprintf('select=%d\n',select);
    %fprintf('n_actions=%d %d %d %d\n',n_actions.');
    %fprintf('n_actions=%d %d %d\n',n_actions.');
    %fprintf('n_actions=%d %d\n',n_actions.');
    temp=-100;
    
    if (select==1) || (((length(n_actions)==4)) && (q(cs,n_actions(1))==0) && (q(cs,n_actions(2))==0) && (q(cs,n_actions(3))==0) && (q(cs,n_actions(4))==0)) || (((length(n_actions)==3)) && (q(cs,n_actions(1))==0) && (q(cs,n_actions(2))==0) && (q(cs,n_actions(3))==0)) || (((length(n_actions)==2)) && (q(cs,n_actions(1))==0) && (q(cs,n_actions(2))==0));
        ns = n_actions(randi([1 length(n_actions)],1,1));
        %fprintf('1st action executed\n');
    else
    
    
    temp11=0;
        for g=1:length(n_actions)
            if q(cs,n_actions(g))>temp
                temp=q(cs,n_actions(g));
                temp11=n_actions(g);
            end
            %fprintf('q(cs,n_actions(g)=%d\n',q(cs,n_actions(g)));
        end
        %fprintf('2nd action executed\n');
        ns=temp11;
    end
    
    
                
    %fprintf('temp=%d',temp);        
    %fprintf('ns=%d',ns);
    
    for x=1:n
        for y=1:n
            %fprintf('x=%d\n)',x);
            %fprintf('y=%d\n',y);
            if ns==smatrix(x,y)
                temp1=x;
                temp2=y;
                %fprintf('temp1=%d\n)',x);
                %fprintf('temp2=%d\n',y);
            end
        end
    end
    
    ns_row=temp1;
    ns_col=temp2;
    
       
            % find all the possible actions for the selected state
            if (ns_row==1)&&(ns_col==1)
                n_actions =[ns_col+1 n+1];
            end
    if (ns_row==1)&&(ns_col>1 && ns_col<n)
        n_actions =[ns_col-1 ns_col+1 ns_col+n];
        
    end
    if (ns_row==1)&&(ns_col==n)
        n_actions =[ns_col-1 2*n];
        
    end
    if (ns_row>1 && ns_row<n) && (ns_col==n)
        n_actions =[n*(ns_row-1) (n*ns_row-1) n*(ns_row+1)];
        
    end
    if (ns_row==n) && (ns_col==n)
        n_actions =[n*(n-1) (n*n-1)];
        
    end
    if (ns_row==n) && (ns_col<n && ns_col>1)
        n_actions =[n*(ns_row-1)+ns_col+1 n*(ns_row-1)+ns_col-1 n*(ns_row-1)+ns_col-n];
        
    end
    if (ns_row==n) && (ns_col==1)
        n_actions =[n*(ns_row-1)+2 n*(ns_row-1)+1-n];
        
    end
    if (ns_row<n && ns_row>1) && (ns_col==1)
        n_actions =[n*(ns_row-1)+2 n*(ns_row-1)+1-n n*(ns_row-1)+1+n];
    end
    if (ns_row<n && ns_row>1) && (ns_col<n && ns_col>1)
        n_actions =[n*(ns_row-1)+ns_col+1 n*(ns_row-1)+ns_col-1 n*(ns_row-1)+ns_col-n n*(ns_row-1)+ns_col+n];
    end
            
            % find the maximum q-value i.e, next state with best action
             max_q = 0;
            for j=1:length(n_actions)
                max_q = max(max_q,q(ns,n_actions(j)));
            end
            
            % Update q-values as per Bellman's equation
            
            
            q(cs,ns)=(1-alpha)*q(cs,ns)+alpha*[reward(ns_row,ns_col)+gamma*max_q];
            cumulative = cumulative + reward(ns_row,ns_col);
            
            if reward(ns_row,ns_col)==-10
                ns=cs;
            end
            
            
            % Check whether the episode has completed i.e Goal has been reached
            if(cs == Goal)
                fprintf('cum_reward=%d\n',cumulative);
                break;
            end
            
            % Set current state as next state
            cs=ns;
            cs_row=ns_row;
            cs_col=ns_col;
    end
    fprintf('Total steps.=%d\n',u);
    step_matrix(i)=u;
    cum_matrix(i)=cumulative;
    epi_matrix(i)=i;
end

figure
plot(epi_matrix,step_matrix,'g');
title('No. of steps vs. Episode no.');
xlabel('Episode no.');
ylabel('No. of Steps');

figure
plot(epi_matrix,cum_matrix,'g');
title('Cumulative reward vs. Episode no.');
xlabel('Episode no.');
ylabel('Cumulative reward');

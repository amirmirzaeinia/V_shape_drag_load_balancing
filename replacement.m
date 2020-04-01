% Fixed size V-shape swarm replacement
% it generates the energy consumption graphs as well as output positioing
% matrix to be used in flying animation module

%%
% prompt = {'Enter flock size: odd number' , 'Enter S', 'a' , 'b' , 'm(kg)' , 'v(mps)' , 'Dist(m)' , 'Init Energy(w)'};
% dlgtitle = 'Input';
% dims = [1 35];
% definput = {'9', '.055' , '0.432' , '.55' , '1.1' , '15', '10' , '54000'};
% answer = inputdlg(prompt,dlgtitle,dims,definput);
% 
% N = str2num(answer{1});
% s = str2num(answer{2});
% a = str2num(answer{3});
% b = str2num(answer{4});
% m = str2num(answer{5});
% v = str2num(answer{6});
% dist = str2num(answer{7});
% init_energy= str2num(answer{8});

%%
% definition of input values
%N is number of agents
%s is the distance between the wing tips
%a is vortex length in the wing tip
%b is the semi wing span
%m is the mass
%v is the speed of the swarm
%dis is the distance the swarm can fly
%init_energy is the amount of energy that each agent has at the begining
%these values can be given as inpput or preset default value.
%%

%%

% *******************************************************************************************
% *******************************************************************************************
% *******************************************************************************************
% the main ffunction will call the whole process
% you may look at this as a black box, just enter the replacement and you will get the
% output based on the defult values
% if you wanna change the defaul values you may pass them to the main
% function and comment out the default value. 

%there are two output matrix as sorted energy and the agent ids
% Following is a sample output 
% Notice that each row of the energy matrix is sorted and their IDs are
% available in ID matrix
% the first row of the energy matrix shows the highest energy in each step
% Another point is that replacement happens once the leader agent get to the half of its
% energy 
% *******************************************************************************************
% *******************************************************************************************
% *******************************************************************************************
%%

%%
% sample output
% Returned_energy_matrix =
% 
%    1.0e+04 *
% 
%     5.4000    3.4304    2.1560    1.2878    0.6439
%     5.4000    3.4304    2.1560    1.2878    0.6439
%     5.4000    3.4162    2.0948    1.0826    0.5960
%     5.4000    3.4162    2.0948    1.0826    0.5960
%     5.4000    3.3788    1.9298    1.0780    0.5765
%     5.4000    3.3788    1.9298    1.0780    0.5765
%     5.4000    3.2779    1.7152    0.9232    0.4501
%     5.4000    3.2779    1.7152    0.9232    0.4501
%     5.4000    2.7000    1.4488    0.6624    0.1927
%     5.4000    2.7000    1.4488    0.6624    0.1927
% 
% 
% id =
% 
%      1     5     7     4     4
%      2     6     2    10     1
%      3     4     4     2     2
%      4     7     9     3     5
%      5     3     6     7    10
%      6     8     1     8     6
%      7     2     5     9     3
%      8     9    10     5     9
%      9     1     3     6     7
%     10    10     8     1     8
    %%
clc
clear all
close all

main

function [Returned_energy_matrix id] = main ()
         N = 14;
        s = -.0825;
        a = 0.432;
        b = .75;
        m = 3.8;
        v = 18;
        dist =  10;
        init_energy = 54000;
        
        d_11 = (2*m^2*9.8^2)/(pi*1.225*b^2*v^2);
        [Returned_energy_matrix id] = flock_size_drag_fun(a, s , b , N , d_11 , v , dist , init_energy)


end
%%




%%
function   [Returned_energy_matrix id] = flock_size_drag_fun(a,s,b,N,d_11 , v , dist , init_energy)


%          title('Energy consuming after replacement')
         
            [ d_array] = tmp_drag_fun (a , s , b , N , d_11 , v , dist , init_energy);   
            
            i = 1;
            E_remind = ones(1 , N) * init_energy;
            Energy_total = ones (1 , N) * .5 * init_energy;
            
            non_replacement_dist = .99 * init_energy/(d_array(1));
            non_replacement_time=  non_replacement_dist / v;
            non_replacement_energy = E_remind;
            id (:,i) = 1:N;
            Returned_energy_matrix(:,i) = init_energy.* ones(N,1);
            while( 1)  
                distance(i) =   .5 * E_remind(1) / d_array(1);
                time(i) = distance(i) / v ;
                
                Energy_total  = distance(i) .* d_array ;
                E_remind =  E_remind - Energy_total;
                non_replacement_energy = non_replacement_energy-Energy_total;
                if ( min ( E_remind ) < .01 * init_energy)
                         break
                end
                   
                figure 
                plot(E_remind)
                hold on    
                
                
                 title('Energy remind')
                tmp_E_remind = E_remind;
                tmporar_E_remind =E_remind;
                
                [max_reminder index_lead] = max(tmp_E_remind);
                tmp_E_remind(index_lead) = 0;
                
                E_remind (index_lead) = E_remind(1);
                E_remind (1) = max_reminder;
                
                [max_reminder index_tail] = max(tmp_E_remind);
                E_remind (index_tail) = E_remind(end);
                E_remind (end) = max_reminder;
                
                

                plot ( E_remind , '-*')
                hold off
                legend('Reminded Energy Before replacement','Reminded Energy after replacement')
                
                tmp_lead_id = index_lead;
                tmp_tail_id = index_tail;
                
                
                i = i + 1;
                
                [tmp_sorted_E_remind,index] = sort(tmporar_E_remind , 'descend');
                id(:,i) = id( index,i-1);
                Returned_energy_matrix(:,i) = tmp_sorted_E_remind;
            
                
            end
            
            
            dis_aggregate=0;
            agg =0;
            for i = 1 : length(distance)
                agg = agg + distance(i);
                dis_aggregate(i) = agg; 
            end
            
            figure
            plot(1:length(dis_aggregate),non_replacement_dist*ones(length(dis_aggregate)))
            hold on
            bar(dis_aggregate)
            title('Total Distance in each replacement')
            
            agg_time=0;
            time_aggregate=0;
            for i = 1:length(time)
                agg_time = agg_time + time(i);
                time_aggregate(i) = agg_time; 
            end
            figure
            plot(1 : length(time_aggregate),non_replacement_time*ones(length(dis_aggregate)))
            hold on
            bar(time_aggregate)
             title('Time')

       
end
%%


%%
function [  d_array]= tmp_drag_fun(a , s , b , N , d_11, v , dist , init_power)

    for i = 1:N
        d_array(i) = drag_fun(a,s,b,i,N,d_11) ; 
    end

    
end
%%

function d_total = drag_fun(a,s,b,j,N,d_11)
    drag_val = 1;

    for i = 1 : N
            if i ~= j
                denum = abs(i-j)*(2*b+s);
                drag_val = drag_val*(1-(2*a/denum).^2);

            end
        
    end
    
    d_total = log(drag_val) * (2 * d_11 )/ (pi ^ 2) +d_11;
    

end
%%


clc
clear all
close all



prompt = {'Enter flock size: odd number','Enter S','a','b','m(kg)','v(mps)','dist(m)' , 'Init power(w)'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'14', '.23' , '0.55' , '1.4' , '1.3' , '14', '54000' , '15'};
answer = inputdlg(prompt,dlgtitle,dims,definput);

N = str2num(answer{1});
s = str2num(answer{2});
a = str2num(answer{3});
b = str2num(answer{4});
m = str2num(answer{5});
v = str2num(answer{6});
dist = str2num(answer{7});
init_power= str2num(answer{8});



d_11 = (2*m^2*9.8^2)/(pi*1.225*b^2*v^2);

flock_size_drag_fun(a,s,b,N,d_11,v,dist,init_power)


function  flock_size_drag_fun(a,s,b,N,d_11 , v , dist , init_power)


         for k = 1:N
                [ d_array] = tmp_drag_fun(a,s,b,N,d_11,v,dist,init_power);   

              x=1:N;
            plot(x,d_array,'-o')
            hold on
            title('Individual UAV drag reduction')
                xlabel('Flock size') 
                ylabel('Drag value')

        end
        figure
         hold on

        for k = 1:N
                [ d_array] = tmp_drag_fun(a,s,b,N,d_11,v,dist,init_power);   

              x=1:N;
              p_total = d_array .* v ;
            plot(x,p_total,'-o')
            title('Consumed power')
                xlabel('Flock size') 
                ylabel('power')

        end
        figure
         hold on
         
         for k = 1:N
                [ d_array] = tmp_drag_fun(a,s,b,k,d_11,v,dist,init_power);   

              x=1:k;
              p_total = d_array .* v ;
              p_remind = init_power - p_total;
            plot(x,p_remind,'-o')
            title('Remined power')
                xlabel('Flock size') 
                ylabel('Power')

         end
        
         
          figure
         hold on
         
         for k = 1:N
                [ d_array] = tmp_drag_fun(a,s,b,k,d_11,v,dist,init_power);   

              x=1:k;
              Energy_consumed = d_array * dist;
              
            plot(x,Energy_consumed,'-o')
            title('Energy consumed')
                xlabel('Flock size') 
                ylabel('Energy')

         end
         

end






function [ d_array]= tmp_drag_fun(a,s,b,N,d_11, v , dist,init_power)

    for i = 1:N
        d_array(i) = drag_fun(a,s,b,i,N,d_11) ; 
%          d_array(i) =  d_array(i) + i * d_11;
    end

    
end


function d_total = drag_fun(a,s,b,j,N,d_11)
    drag_val = 1;

    for i = 1:N
            if i ~= j
                denum = abs(i-j)*(b+s);
                drag_val = drag_val*(1-(2*a/denum).^2);

            end
        
    end
    
    d_total = (2*d_11/ (pi^2)) * log(drag_val) + d_11;
    

end



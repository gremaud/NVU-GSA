from single_evaluation_plotting import single_eval_with_plot

values='Robin_Pre2'
model='pre'

list_to_plot=range(1,59+1)



for i in range(len(list_to_plot)):
    for j in range(i+1,len(list_to_plot)):
        removed=[list_to_plot[i],list_to_plot[j]]
        single_eval_with_plot(values,model,removed)
        print(removed)
        
list_to_plot=range(59+1)
for i in range(len(list_to_plot)):
        removed=list_to_plot[i]
        single_eval_with_plot(values,model,removed)
        print(removed)
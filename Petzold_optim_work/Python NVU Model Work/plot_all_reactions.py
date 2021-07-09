from single_evaluation_plotting import single_eval_with_plot

values='Robin_Pre2'
model='pre'

#list_to_plot=range(59+1)
list_to_plot=[[4,52],[4,11],[11,52],[4,23],[23,52],[11,23],[4,13],
              [11,13],[13,52],[13,23],[4,25],[8,11],[11,25]]

for i in range(len(list_to_plot)):
    removed=list_to_plot[i]
    single_eval_with_plot(values,model,removed)
    print(removed)
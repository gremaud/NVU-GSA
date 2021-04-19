from single_evaluation_plotting import single_eval_with_plot

values='Robin_Pre2'
model='pre'

for i in range(70+1):
    single_eval_with_plot(values,model,i)
    print(i)
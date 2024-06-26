import matplotlib.pyplot as plt
import numpy as np

# Data preparation
model_names = ["MIP", "MIP + L2", "GERTHEIS"]
models = [
    [8.89275014813345, 3.410710184518646, 2.4088934786331357, 7.027182824458129, 7.218430707368456, 0.0, 0.0, 0.0, 0.0, 0.0],
    [1.243833773093219, 1.3559798550230329, 1.2450843738727446, 1.1659868534919755, 1.4268669554062245, 0.0, 0.0, 0.0, 0.0, 0.0],
    [2.34200029457317, 2.285715594990576, 1.7210310709646381, 0.9930586606578824, 1.557881949512008, 0.678743818136545, 0.8101792666353592, 0.8438125762959849, 0.28278941773652727, 0.758967810916494]
]

# Plot candle plots for squared errors of each model for each predictor
num_predictors = len(models[0])
num_subplots = num_predictors if num_predictors <= 6 else 6

plt.figure(figsize=(12, 8))
for i in range(num_subplots):
    plt.subplot(2, 3, i+1)
    for j, model in enumerate(models):
        color = plt.cm.tab10(j)  # Choose a unique color for each model
        if(j<5):

            y_values = [0.01, model[i]/2,model[i],model[i]*1.2]  # Wrap the single value in a list
        else:
            y_values = [model[i]]
        x_values = [j+1]  # Position for this box
        plt.boxplot(y_values, positions=x_values, showmeans=False, meanline=True, patch_artist=True, boxprops=dict(facecolor=color), medianprops=dict(linewidth=3, color="black", linestyle='solid'))

    plt.xticks(np.arange(1, len(model_names)+1), model_names)
    plt.ylim(0, None)  # Ensure y-axis starts from 0
    plt.xlabel("Model")
    plt.ylabel(f"Predictor {'6-10' if i+1 >5 else i+1  } SE")
    plt.title(f"Predictor {'6-10' if i+1 >5 else i+1  } SE across Models")
    plt.grid(True)

plt.tight_layout()
plt.show()

import pandas as pd
import json

# Provided data
# Load data from JSON file
sim = "sim3"
type = "lambda"
ols= 0
# Load data from JSON file
with open(f'{sim}/{type}/results.json', 'r') as f:
    data = json.load(f)
# Initialize an empty DataFrame
summary_df = pd.DataFrame(columns=["Model","#Basis","Lambda",  "ISE_Coefficients", "RMSE_Predictions", "MSE_Predictions"])#,"IQR","R_squared"])

# Iterate over the data and populate the DataFrame
for group, run in (data.items()):
    # loop between each run
    for i, run_data in enumerate(run):
        # if(ols == 0 and i == 1):
        #     break
        summary_df = summary_df._append({
            "#Basis": 10,
            "Lambda": group,
            "Model":  "Fold" + str(i+1),
            "ISE_Coefficients": run_data["ISE_Coefficients"],
            "RMSE_Predictions": run_data["RMSE_Predictions"],
            "MSE_Predictions": run_data["MSE_Predictions"],
            # "IQR": run_data["IQR"],
            #  "R_squared" : run_data["R_squared"]

        }, ignore_index=True)
   
   
# Sort the DataFrame by the "Lambda" column
summary_df["Lambda"] = pd.to_numeric(summary_df["Lambda"])
summary_df = summary_df.sort_values(by="Lambda")
       

# Display the summary table
print(summary_df)
save_path = f'{sim}/{type}\summary_table.tex'
summary_df.to_latex(save_path, index=False)
summary_df.to_html(f'{sim}/{type}\summary_table.html', index=False)
summary_df.to_excel(f'{sim}/{type}\summary_table.xlsx', index=False)
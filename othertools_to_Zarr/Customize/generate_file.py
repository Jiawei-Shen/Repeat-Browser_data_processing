import pandas as pd

def main():
    # Create a sample DataFrame
    data = {
        'TE': ['MER41B', 'LTR13', 'LTR12'],
        'uni_value': [25, 30, 35],
        'all_value': [25, 30, 35],
    }

    df = pd.DataFrame(data)

    # Write the DataFrame to a CSV file
    df.to_csv('./sample_files/heatmap.csv', index=False)  # Set index=False to exclude the index column


if __name__ == "__main__":
    main()
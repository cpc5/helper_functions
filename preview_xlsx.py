import pandas as pd
from pathlib import Path

#### preview excel files via pandas
#### Carlos Perez-Cervantes. Moskowitz lab

def preview_excel(file_path, rows=5, cols=None, sheet_idx=1):
    """
    Preview an Excel file in the terminal.
    
    Args:
        file_path (str): Path to the Excel file
        rows (int): Number of rows to display (default: 5)
        cols (int): Number of columns to display (if None, shows all columns)
        sheet_idx (int): Sheet index to display (default: 1, first sheet)
    """
    try:
        # Check if file exists
        if not Path(file_path).exists():
            print(f"Error: File '{file_path}' not found.")
            return
            
        pd.set_option('display.max_columns', cols)  # Set maximum columns to display
        
        # Get list of sheets
        xl = pd.ExcelFile(file_path)
        num_sheets = len(xl.sheet_names)
        
        # Validate sheet index
        if sheet_idx < 1 or sheet_idx > num_sheets:
            print(f"Error: Sheet index must be between 1 and {num_sheets}")
            print(f"Available sheets:")
            for i, name in enumerate(xl.sheet_names, 1):
                print(f"{i}. {name}")
            return
            
        # Read specified sheet (converting from 1-based to 0-based index)
        df = pd.read_excel(file_path, sheet_name=sheet_idx-1)
        
        print(f"\nFile: {file_path}")
        print(f"Sheet {sheet_idx} of {num_sheets}: '{xl.sheet_names[sheet_idx-1]}'")
        print(f"Shape: {df.shape[0]} rows × {df.shape[1]} columns")
        print("\nColumn names:")
        for i, col in enumerate(df.columns, 1):
            print(f"{i}. {col}")
        print(f"\nFirst {rows} rows:")
        if cols:
            print(df.iloc[:rows, :cols])
        else:
            print(df.head(rows))
            
    except Exception as e:
        print(f"Error reading Excel file: {str(e)}")

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Preview Excel files in terminal')
    parser.add_argument('file', help='Path to Excel file')
    parser.add_argument('-r', '--rows', type=int, default=5, help='Number of rows to display')
    parser.add_argument('-c', '--cols', type=int, help='Number of columns to display')
    parser.add_argument('-s', '--sheet', type=int, default=1, help='Sheet index to display (1-based)')
    
    args = parser.parse_args()
    preview_excel(args.file, args.rows, args.cols, args.sheet)

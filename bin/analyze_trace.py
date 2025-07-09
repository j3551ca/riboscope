#!/usr/bin/env python3
#%%
import pandas as pd
import sys
import os
from datetime import datetime
import argparse
import csv
#%%
dir = ""
def load_trace_csv(trace_file):
    """Load and validate trace data
    """
    try:
        df = pd.read_csv(trace_file, sep='abcxyz', header=0)
        print(f"Loaded {len(df)} tasks from {trace_file}")
        return df
    except Exception as e:
        print(f"Error loading trace file: {e}")
        return None
    
def load_trace_data(trace_file):
    """Load and validate trace data
    """
    try:
        records = []
        current_record = []
        header = None
        found_header = False

        with open(trace_file, "r") as f:
            for line in f:
                # new record starts when the line does NOT begin with whitespace - only diff between new record lines and other 
                if line.strip() and not line.startswith((' ', '\t', '"', ">", "'")):
                    if not found_header:
                # assume first non-indented line is the header
                        header = line.strip().split("abcxyz")
                        found_header = True
                        continue
                    elif current_record: #prev record exists
                            records.append(''.join(current_record)) # smush script lines together with current record kept in records because field sep is unique
                            current_record = [] # next record/ actual real line
                    current_record.append(line) # start record
                elif line.strip().startswith("#") or line.strip()=="":
                    continue # skip commented lines and empty new lines
                else:
                    current_record.append("SCRIPTNEWLINE" + line) # script field lines that do start with above characters
            if current_record:
                records.append(''.join(current_record))

        parsed = []
        for rec in records:
            # remove newlines, normalize spacing
            flattened = rec.replace('\n', ' ').replace('\r', ' ')
            fields = flattened.split("abcxyz")
            if len(fields) != len(header):
                print(f"Skipping malformed line - expected {len(header)}, but got {len(fields)}")
                continue
            entry = { k.strip(): v.strip() for k,v in zip(header, fields)}
            
            parsed.append(entry)
        
        df = pd.DataFrame(parsed)
        return df
                
    except Exception as e:
        print(f"Error loading trace file: {e}")
        return None
    
def sanitize_sample_ids(df):
    """Eliminate non-sample ID characters from tag
    """
    df["tag"] = df["tag"].str.split().str[0]
    return df

#%%
trace_file=os.path.join(dir,"execution_trace.csv")

df = load_trace_data(trace_file)
trace_df = sanitize_sample_ids(df)


#%%

# parsing script field:
def parse_bash_commands(command_string: str) -> dict:
    """
    Parses a bash command string into its tool, parameters, and positional files.
    Args:
        command_string: The raw command string from Nextflow.
    Returns:
        A dictionary with 'tool', 'params', and 'files'.
    """
    if not isinstance(command_string, str):
        return {
            "tool": None,
            "params": {},
            "files": []
        }
    try:
        # split into tokens with shlex - for case when quoted subcommands and split() would separate them
        tokens = shlex.split(command_string)
    except ValueError:
        return { # if it doesn't parse properly, store it as-is
            "tool": command_string,
            "params": {},
            "files": []
        }
    
    if not tokens:
        return {"tool": None, 
                "params": {}, 
                "files": []}
    
    print(tokens)
    params = {}
    remaining_tokens = []
    i = 0
    while i < len(tokens):
        token = tokens[i]
        if token.startswith("-"):
            if "=" in token:
                flag, value = token.split("=",1)
                params[flag] = value
                i +=1
                continue
            flag = token
            # is next element value or is boolean flag
            if (i + 1) < len(tokens) and not tokens[i + 1].startswith("-"):
                params[flag]=tokens[i + 1]
                i +=2 # skip the flag and value tokens to next 
            else:
                params[flag] = True
                i +=1
        else:
            remaining_tokens.append(token)
            i +=1
    
    print(params)
    print(remaining_tokens)

    # assuming anything with a . is a file - could be fishy
    files = []
    files = [ token for token in remaining_tokens if "." in token]
    files = list(set(files))
    print(files)
    commands = [ token for token in remaining_tokens if token not in files and token.isalpha() ]

    print(commands)

    return({"commands": commands,
           "parameters":params,
           "files": files})

#%%
# 600
parse_bash_commands(trace_df.iloc[123]["script"])


#%%

import shlex

def parse_bash_commands(command_string: str) -> dict:
    """
    Parses a bash command string into its tool, parameters, and positional files.

    Args:
        command_string: The raw command string from Nextflow.

    Returns:
        A dictionary with 'tool', 'params', and 'files'.
    """
    if not isinstance(command_string, str):
        return {
            "tool": None,
            "params": {},
            "files": []
        }
    try:
        # 1. Split the command using shlex for robustness
        tokens = shlex.split(command_string)
    except ValueError:
        return { # Handle cases of malformed strings, e.g., unmatched quotes
            "tool": command_string,
            "params": {},
            "files": []
        }

    if not tokens:
        return {
            "tool": None,
            "params": {},
            "files": []
        }

    # 2. Identify the base tool command
    tool_parts = []
    token_index = 0
    while token_index < len(tokens) and not tokens[token_index].startswith('-'):
        tool_parts.append(tokens[token_index])
        token_index += 1
    
    tool = " ".join(tool_parts)

    # 3. Parse parameters and positional files from the rest of the tokens
    params = {}
    files = []
    remaining_tokens = tokens[token_index:]
    i = 0
    while i < len(remaining_tokens):
        token = remaining_tokens[i]
        if token.startswith('-'):
            # Handle '--key=value' format
            if '=' in token:
                key, value = token.split('=', 1)
                params[key] = value
                i += 1
                continue

            key = token
            # Check if the next token is a value or another flag
            if (i + 1) < len(remaining_tokens) and not remaining_tokens[i+1].startswith('-'):
                params[key] = remaining_tokens[i+1]
                i += 2 # Move past both the key and the value
            else:
                # It's a boolean flag (e.g., --verbose)
                params[key] = True
                i += 1
        else:
            # It's a positional argument/file
            files.append(token)
            i += 1
            
    return {
        "tool": tool,
        "params": params,
        "files": files
    }
#%%
# --- Apply the function to the DataFrame column ---
# The result of .apply() will be a Series of dictionaries.
parsed_script = trace_df['script'].apply(parse_bash_commands)

#%%
# You can convert this Series of dicts into a new DataFrame
final = pd.json_normalize(parsed_script)

# Combine with the original df for a full picture
final_df = pd.concat([trace_df, final], axis=1)

#%%
def main(args):
    df = load_trace_data(args.trace_file)
    df.to_csv(os.path.join(dir, "analyzed_trace.csv"), index=False)
    if df is None:
        sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Analyze Nextflow execution trace')
    parser.add_argument('--trace_file', required=True, help='Path to execution trace CSV file')
    args = parser.parse_args()
    main(args)

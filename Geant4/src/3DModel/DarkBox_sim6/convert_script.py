import os

def force_clean_for_ubuntu(file_path):
    # Read text directly
    with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()

    if not lines:
        return False

    # Force the first line to be exactly 'solid ASCII' without spaces or carriage returns
    first_line = lines[0].strip()
    
    # Check if it contains 'solid' regardless of case
    if 'solid' in first_line.lower():
        # Clean lines and remove carriage returns (\r) from all lines
        cleaned_lines = []
        
        # Enforce exact clean string initialization for line 1
        cleaned_lines.append("solid ASCII\n")
        
        for line in lines[1:]:
            cleaned_lines.append(line.replace('\r', '').replace('\n', '\n'))

        # Write out using pure raw ascii encoding
        print(f"🛠️ Re-encoding and aligning stream for: {os.path.basename(file_path)}")
        with open(file_path, 'w', encoding='ascii') as f:
            f.writelines(cleaned_lines)
        return True

    return False

def process_directory(directory="."):
    print(f"🔄 Re-formatting directory: '{os.path.abspath(directory)}'\n")
    count = 0
    for filename in os.listdir(directory):
        if filename.lower().endswith('.stl'):
            if force_clean_for_ubuntu(os.path.join(directory, filename)):
                count += 1
    print(f"\n✅ Cleaned and stream-aligned {count} STL files.")

if __name__ == "__main__":
    process_directory()

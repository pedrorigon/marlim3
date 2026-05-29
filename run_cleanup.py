import os
import json
import glob

BASE = r'c:\Users\cvpo\GitHub\mrtr\motor'

def remove_keys_from_json_file(filepath, keys_to_remove):
    """Remove specified top-level keys from a JSON file."""
    with open(filepath, 'r', encoding='utf-8') as f:
        data = json.load(f)
    
    changed = False
    for key in keys_to_remove:
        if key in data:
            del data[key]
            changed = True
    
    if changed:
        with open(filepath, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=2, ensure_ascii=False)
            f.write('\n')
        print(f"  Updated: {filepath}")
    return changed

# MVM JSON files only (other files already handled)
print("=== MVM JSON files ===")
mvm_files = glob.glob(os.path.join(BASE, 'mvm', 'mr3', '**', '*.json'), recursive=True)
print(f"  Found {len(mvm_files)} mvm JSON files")
count = 0
errors = 0
for f in mvm_files:
    try:
        if remove_keys_from_json_file(f, ['jsonVersion', 'versaoJson', 'versaoJSON', 'versao', 'version']):
            count += 1
    except Exception as e:
        print(f"  ERROR in {f}: {e}")
        errors += 1
print(f"  Updated {count} mvm files, {errors} errors")
print("Done!")

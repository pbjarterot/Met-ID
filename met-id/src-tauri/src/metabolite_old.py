import sys
import json
import struct
from rdkit import Chem

def substruct_matches(smiles_strings, smarts_pattern):
    smarts_mol = Chem.MolFromSmarts(smarts_pattern)
    if smarts_mol is None:
        return "Invalid SMARTS pattern"
    
    results = []
    for smiles in smiles_strings:
        input_mol = Chem.MolFromSmiles(smiles)
        if input_mol is None:
            results.append(0)
        else:
            matches = len(input_mol.GetSubstructMatches(smarts_mol))
            results.append(matches)
    return results

def main():
    # Read the length of the incoming data
    length_data = sys.stdin.buffer.read(4)
    length = struct.unpack('I', length_data)[0]

    # Read the actual data
    input_data = sys.stdin.buffer.read(length)
    
    # Deserialize the JSON to a list of strings
    data = json.loads(input_data)
    #print(f"Child received data length: {len(data)}")

    # Process the data to integers
    result = substruct_matches(data[1:], data[0])  # Example: convert strings to their lengths

    # Serialize the result to JSON
    result_json = json.dumps(result)
    result_bytes = result_json.encode()

    # Send the length of the result followed by the actual result
    sys.stdout.buffer.write(struct.pack('I', len(result_bytes)))
    sys.stdout.buffer.write(result_bytes)
    sys.stdout.buffer.flush()

if __name__ == "__main__":
    main()
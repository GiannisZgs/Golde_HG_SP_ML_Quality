import scipy.io
import json
import numpy as np
import os


def convert_nested_matlab_to_json(mat_file_path, json_file_path=None):
    """
    Convert a MATLAB .mat file with nested structures to a JSON file.
    
    Parameters:
    -----------
    mat_file_path : str
        Path to the input .mat file
    json_file_path : str, optional
        Path for the output JSON file. If None, it will use the same name as the input file
        but with a .json extension
    
    Returns:
    --------
    str
        Path to the created JSON file
    """
    
    # If no output path specified, create one based on input path
    if json_file_path is None:
        json_file_path = os.path.splitext(mat_file_path)[0] + '.json'
    
    # Load the MATLAB file
    print(f"Loading MATLAB file: {mat_file_path}")
    mat_data = scipy.io.loadmat(mat_file_path, squeeze_me=True, struct_as_record=False)
    
    # Clean up MATLAB-specific entries
    for key in list(mat_data.keys()):
        if key.startswith('__'):
            del mat_data[key]
    
    print("Converting MATLAB data to JSON-compatible format...")
    
    # Custom JSON encoder to handle numpy arrays and MATLAB objects
    class MatlabStructEncoder(json.JSONEncoder):
        def default(self, obj):
            # Handle numpy arrays
            if isinstance(obj, np.ndarray):
                if obj.ndim == 0:
                    return obj.item()
                return obj.tolist()
            
            # Handle MATLAB objects (loaded as scipy.io.matlab.mio5_params.mat_struct)
            elif isinstance(obj, scipy.io.matlab.mio5_params.mat_struct):
                return {attr: getattr(obj, attr) for attr in dir(obj) if not attr.startswith('__')}
            
            # Handle numpy data types
            elif isinstance(obj, (np.int_, np.intc, np.intp, np.int8, np.int16, np.int32, np.int64,
                                 np.uint8, np.uint16, np.uint32, np.uint64)):
                return int(obj)
            elif isinstance(obj, (np.float_, np.float16, np.float32, np.float64)):
                return float(obj)
            elif isinstance(obj, (np.bool_)):
                return bool(obj)
            elif isinstance(obj, complex):
                return str(obj)
            elif isinstance(obj, bytes):
                print(f"Encountered bytes object: {obj}")
                try:
                    return obj.decode('utf-8')  # Attempt to decode as UTF-8
                except UnicodeDecodeError:
                    return obj.hex()  # Fallback: represent as a hex string

            # Let the base class default method handle the rest
            return json.JSONEncoder.default(self, obj)
    
    # Save to JSON
    try:
        with open(json_file_path, 'w') as json_file:
            json.dump(mat_data, json_file, cls=MatlabStructEncoder, indent=2)
        print(f"Successfully converted to JSON: {json_file_path}")
        return json_file_path
    except Exception as e:
        print(f"Error saving JSON file: {e}")
        raise

if __name__ == "__main__":
   
    input_file = "/home/giannis/Documents/ECG HG paper/results_data/metrics_deviation_from_noise.mat"
    output_file = "/home/giannis/Documents/ECG HG paper/results_data/metrics_deviation_from_noise.json"
    
    convert_nested_matlab_to_json(input_file, output_file)
from flask import Flask, request, jsonify
import pandas as pd
import pickle
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
from flask_cors import CORS  # For handling Cross-Origin Resource Sharing
from sklearnex import patch_sklearn
patch_sklearn()


# Load the pre-trained BBBP model
model = pickle.load(open('hepg2_model.pkl', 'rb'))

# Featurizer class for molecular fingerprint generation (as defined earlier)
class Featurizer:
    def __init__(self, y_column=None, smiles_col='Drug', **kwargs):
        self.y_column = y_column
        self.smiles_col = smiles_col
        self.__dict__.update(kwargs)

    def __call__(self, df):
        raise NotImplementedError()

class ECFPFeaturizer(Featurizer):
    def __init__(self, y_column=None, radius=2, length=1024, **kwargs):
        self.radius = radius
        self.length = length
        super().__init__(y_column, **kwargs)

    def __call__(self, df):
        fingerprints = []
        labels = []
        for i, row in df.iterrows():
            smiles = row[self.smiles_col]
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError(f"Invalid SMILES string: {smiles}")
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, self.radius, nBits=self.length)
            fingerprints.append(np.array(fp))  # Convert to NumPy array

        fingerprints = np.array(fingerprints)  # Create NumPy array of fingerprints
        return fingerprints, None

# Initialize Flask application
app = Flask(__name__)
CORS(app, resources={r"/*": {"origins": "*"}})  # Enable CORS for all routes

# Initialize the featurizer (same as used in model training)
featurizer = ECFPFeaturizer(smiles_col='Drug')

# Route to check API status
@app.route('/', methods=['GET'])
def get_data():
    data = {
        "message": "HEPG2 Model API is Running"
    }
    return jsonify(data)

# Define the route for making predictions
@app.route('/predict', methods=['POST'])
def predict():
    try:
        # Get the SMILES string from the POST request
        data = request.get_json()
        print(f"Received data: {data}")  # Debugging line to see incoming request
        smiles_string = data.get('smiles')

        if not smiles_string:
            return jsonify({'error': 'No SMILES string provided'}), 400

        # Create a DataFrame with the SMILES string
        df = pd.DataFrame({'Drug': [smiles_string]})
        X_new, _ = featurizer(df)

        # Ensure that the input is reshaped as needed
        if X_new.shape[1] == 0:
            return jsonify({"error": "Input features are empty. Please provide a valid SMILES string."}), 400

        X_new = X_new.reshape(1, -1)  # Reshape to be 2D

        # Model prediction
        prediction_class = model.predict(X_new)[0]
        prediction_proba = model.predict_proba(X_new)[:, 1][0]

        # Return the prediction as a JSON response
        return jsonify({
            'Predicted Class': int(prediction_class),
            'Predicted Probability for Class 1': float(prediction_proba)
        })

    except Exception as e:
        # Handle any exceptions and return the error message
        return jsonify({'error': str(e)}), 500

# Run the Flask application
if __name__ == '__main__':
    app.run(debug=True, port=5007)
    # Use host='0.0.0.0' to allow access from other devices on the network
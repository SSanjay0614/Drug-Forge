from flask import Flask, request, jsonify
from flask_cors import CORS
import pickle
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

app = Flask(__name__)
CORS(app)  # Enable CORS for frontend communication

# Load the trained model
with open('solubility_model.pkl', 'rb') as f:
    model = pickle.load(f)

# Function to featurize SMILES string
def featurize_smiles(smiles, radius=2, length=1024):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fingerprint = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=length)
    return np.array(fingerprint).reshape(1, -1)

@app.route('/', methods=['GET'])
def home():
    return jsonify({'message': 'Welcome to the Solubility Prediction API! Use the /predict endpoint to make predictions.'})

@app.route('/predict', methods=['POST'])
def predict():
    data = request.get_json()
    smiles = data.get('smiles')

    # Featurize the input SMILES string
    features = featurize_smiles(smiles)
    if features is None:
        return jsonify({'error': 'Invalid SMILES string'}), 400

    # Predict using the loaded model
    prediction = model.predict(features)

    return jsonify({'prediction': prediction.tolist()})

if __name__ == '__main__':
    app.run(debug=True)
